#!/usr/bin/env python
"""
scripts/update_zenodo.py -- Automates updating an existing Zenodo deposition with a new release version
using the new InvenioRDM REST API (/api/records).
"""

import argparse
import copy
import datetime
import hashlib
import json
import os
import re
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Iterable, Optional, Union


def parse_pyproject_version(
    filepath: str = "moleditpy/pyproject.toml",
) -> Optional[str]:
    """Extract version from pyproject.toml."""
    if not os.path.exists(filepath):
        return None
    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()
    match = re.search(r'(?m)^version\s*=\s*["\']([^"\']+)["\']', content)
    if match:
        return match.group(1)
    return None


def load_config(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    config = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise ValueError("Zenodo configuration must be a JSON object.")
    for key in ("deposition_id", "sandbox_deposition_id"):
        value = config.get(key)
        if value is not None and (not isinstance(value, str) or not value.isdecimal()):
            raise ValueError(f"{key} must be a numeric record ID string.")
    return config


CONFIG = load_config(Path(__file__).with_name("zenodo-config.json"))
DEPOSITION_ID = CONFIG.get("deposition_id")


def file_md5(path: str) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def plan_file_uploads(
    existing_entries: Union[Iterable[dict[str, Any]], dict[str, dict[str, Any]]],
    paths: list[str],
) -> tuple[list[str], list[str]]:
    """Which local files still need uploading to a draft, and which draft
    entries are stale and must be deleted first.

    existing_entries: the draft's file entries ({"key", "checksum":
    "md5:<hex>", "status"}). A file already committed with the same MD5 is
    skipped; an entry under the same name that differs or never finished
    uploading is stale.
    """
    if isinstance(existing_entries, dict):
        existing_entries = existing_entries.values()
    by_key = {e.get("key"): e for e in existing_entries}
    to_upload, stale = [], []
    for path in paths:
        entry = by_key.get(os.path.basename(path))
        if entry is None:
            to_upload.append(path)
        elif (
            entry.get("status") == "completed"
            and entry.get("checksum") == f"md5:{file_md5(path)}"
        ):
            continue
        else:
            stale.append(entry["key"])
            to_upload.append(path)
    return to_upload, stale


def make_request(
    url: str,
    data: Any = None,
    headers: Optional[dict[str, str]] = None,
    method: str = "GET",
    json_response: bool = True,
) -> Any:
    """Utility function to make HTTP requests using urllib."""
    headers = dict(headers or {})
    is_files = re.search(r"/files(?:/|$)", urllib.parse.urlsplit(url).path)
    headers.setdefault(
        "Accept",
        "application/json" if is_files else "application/vnd.inveniordm.v1+json",
    )

    req_data = None
    if data is not None:
        if isinstance(data, (dict, list)):
            req_data = json.dumps(data).encode("utf-8")
            if "Content-Type" not in headers:
                headers["Content-Type"] = "application/json"
        else:
            req_data = data  # raw bytes

    req = urllib.request.Request(url, data=req_data, headers=headers, method=method)

    try:
        with urllib.request.urlopen(req, timeout=60) as response:
            resp_content = response.read()
            if json_response:
                return json.loads(resp_content.decode("utf-8"))
            return resp_content
    except urllib.error.HTTPError as e:
        err_body = e.read().decode("utf-8", errors="replace")
        req_body_str = ""
        if data is not None and isinstance(data, (dict, list)):
            req_body_str = f"\nRequest Body: {json.dumps(data, indent=2)}"
        raise RuntimeError(
            f"HTTP Error {e.code}: {e.reason}\nURL: {url}\nMethod: {method}{req_body_str}\nResponse Body: {err_body}"
        ) from e
    except urllib.error.URLError as e:
        raise RuntimeError(
            f"URL Error: {e.reason}\nURL: {url}\nMethod: {method}"
        ) from e


def build_draft_update(
    draft: dict[str, Any], version: str, publication_date: str
) -> dict[str, Any]:
    """Preserve editable draft data without copying parent or system fields."""
    metadata = copy.deepcopy(draft.get("metadata", {}))
    if not metadata.get("title") or not metadata.get("creators"):
        raise ValueError("Incomplete RDM draft metadata; refusing to overwrite it.")
    if any("person_or_org" not in c for c in metadata["creators"]):
        raise ValueError("Legacy draft representation received; RDM JSON is required.")
    metadata["version"] = version
    metadata["publication_date"] = publication_date
    # RDM RightsSchema accepts either a vocabulary ID or free text, never both.
    if "rights" in metadata:
        metadata["rights"] = [
            {"id": right["id"]} if right.get("id") else right
            for right in metadata["rights"]
        ]
    payload: dict[str, Any] = {"metadata": metadata, "files": {"enabled": True}}
    access = draft.get("access", {})
    if access:
        payload["access"] = {
            k: copy.deepcopy(access[k])
            for k in ("record", "files", "embargo")
            if k in access
        }
    doi = draft.get("pids", {}).get("doi")
    if doi:
        payload["pids"] = {"doi": copy.deepcopy(doi)}
    if "custom_fields" in draft:
        payload["custom_fields"] = copy.deepcopy(draft["custom_fields"])
    return payload


def check_published_files(record: dict[str, Any], paths: list[str]) -> None:
    """Refuse a duplicate publication even when existing assets differ."""
    entries = record.get("files", {}).get("entries", {})
    for path in paths:
        name = os.path.basename(path)
        entry = entries.get(name, {})
        if entry.get("checksum") != "md5:" + file_md5(path):
            raise ValueError(
                f"Version already published, but asset {name} is missing or differs. "
                "Inspect the published record; no draft was created."
            )


def find_published_version(
    latest: dict[str, Any], version: str, base_url: str, headers: dict[str, str]
) -> Optional[dict[str, Any]]:
    """Find a published version before any mutation, including older backfills."""
    if latest.get("metadata", {}).get("version") == version:
        return latest
    url = f"{base_url}/records/{latest['id']}/versions?size=100"
    seen = set()
    while url:
        if url in seen:
            raise ValueError(
                "Version history pagination repeated; refusing to create a draft."
            )
        seen.add(url)
        page = make_request(url, headers=headers)
        hits = page.get("hits", {}).get("hits")
        if not isinstance(hits, list):
            raise TypeError(
                "Could not verify published version history; no draft created."
            )
        for record in hits:
            if record.get("is_published") is False:
                continue
            if record.get("metadata", {}).get("version") == version:
                published_record: dict[str, Any] = make_request(
                    f"{base_url}/records/{record['id']}", headers=headers
                )
                return published_record
        url = page.get("links", {}).get("next")
    return None


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Update Zenodo Deposition (InvenioRDM API)"
    )
    parser.add_argument("--token", help="Zenodo access token")
    parser.add_argument(
        "--deposition-id",
        default=DEPOSITION_ID,
        help="Zenodo record ID to add a version to (defaults to DEPOSITION_ID)",
    )
    parser.add_argument(
        "--draft-id",
        type=int,
        help="Resume this existing draft without creating another version.",
    )
    parser.add_argument(
        "--version-string",
        help="Version to release (defaults to pyproject.toml version)",
    )
    parser.add_argument(
        "--sandbox",
        action="store_true",
        help="Use Zenodo Sandbox (sandbox.zenodo.org) instead of production",
    )
    parser.add_argument(
        "--publish",
        action="store_true",
        help="Automatically publish the draft deposition after updating",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Perform a dry run (log actions without making API calls)",
    )
    parser.add_argument(
        "files",
        nargs="*",
        help="List of files to upload to the new version deposition",
    )

    args = parser.parse_args()

    # Determine base API URL
    base_url = (
        "https://sandbox.zenodo.org/api" if args.sandbox else "https://zenodo.org/api"
    )

    # Determine token
    token = args.token or os.environ.get("ZENODO_TOKEN")
    if not token and not args.dry_run:
        raise ValueError(
            "Zenodo access token is required. Use --token, --dry-run, or set the ZENODO_TOKEN env variable."
        )

    # Determine version
    version = args.version_string or parse_pyproject_version()
    if not version:
        raise ValueError("Could not determine version. Specify --version-string.")

    # Normalize version: strip leading v/V if present
    if version.lower().startswith("v"):
        version = version[1:]

    deposition_id = args.deposition_id
    if args.sandbox and deposition_id == DEPOSITION_ID:
        deposition_id = CONFIG.get("sandbox_deposition_id")
        if not deposition_id:
            raise ValueError(
                "Pass an explicit Sandbox --deposition-id or configure sandbox_deposition_id."
            )
    if not deposition_id:
        raise ValueError("No deposition ID: pass --deposition-id.")

    print(f"Target Version: {version}")
    print(f"Environment: {'Sandbox' if args.sandbox else 'Production'}")
    print(f"Base API: {base_url}")
    print(f"Deposition ID: {deposition_id}")
    print(f"Files to upload: {args.files}")

    headers = {"Authorization": f"Bearer {token}"} if token else {}

    if args.dry_run:
        print("\n--- DRY RUN: Listing Actions ---")
        if args.draft_id:
            print(f"[DRY-RUN] Would resume existing draft {args.draft_id}.")
        else:
            print("[DRY-RUN] Would check published versions before creating a draft.")
        print("[DRY-RUN] Would fetch draft details using GET on the draft link.")
        if args.files:
            print(
                f"[DRY-RUN] Would register {len(args.files)} files via POST on draft files endpoint."
            )
            for fpath in args.files:
                print(f"[DRY-RUN] Would upload file: {fpath} via PUT content endpoint.")
                print(f"[DRY-RUN] Would commit file: {fpath} via POST commit endpoint.")
        print(f"[DRY-RUN] Would update metadata with PUT request: version={version}")
        if args.publish:
            print("[DRY-RUN] Would publish draft via POST publish endpoint.")
        print("--- DRY RUN COMPLETE ---")
        return

    # Reject invalid local input before creating or modifying a remote draft.
    names = [os.path.basename(path) for path in args.files]
    if len(names) != len(set(names)):
        raise ValueError("Upload filenames must be unique.")
    for path in args.files:
        if not os.path.isfile(path):
            raise FileNotFoundError(path)

    # A fixed record ID may be years behind the current published version.
    parent_url = f"{base_url}/records/{deposition_id}/versions/latest"
    latest = make_request(parent_url, headers=headers)
    published = find_published_version(latest, version, base_url, headers)
    if published is not None:
        check_published_files(published, args.files)
        print(
            f"Version {version} is already published: "
            f"{published.get('links', {}).get('self_html', published.get('id'))}"
        )
        return
    deposition_id = latest["id"]

    # 1. Create a new version of the latest published record.
    print(f"\n[1/6] Creating new version draft for record {deposition_id}...")
    new_version_url = f"{base_url}/records/{deposition_id}/versions"
    if args.draft_id:
        resp = make_request(
            f"{base_url}/records/{args.draft_id}/draft", headers=headers
        )
        parent_id = latest.get("parent", {}).get("id")
        if not parent_id or resp.get("parent", {}).get("id") != parent_id:
            raise ValueError("The selected draft belongs to a different record family.")
    else:
        resp = make_request(new_version_url, data=None, headers=headers, method="POST")

    draft_id = resp.get("id")
    draft_links = resp.get("links", {})
    draft_self_url = draft_links.get("self")
    draft_files_url = draft_links.get("files")
    draft_publish_url = draft_links.get("publish")

    if not draft_id or not draft_self_url:
        raise RuntimeError(f"Could not retrieve draft ID or link from response: {resp}")

    print(f"Draft ID: {draft_id}")
    print(f"Draft URL: {draft_self_url}")
    print(f"Draft Files URL: {draft_files_url}")

    # 2. Fetch the draft details to inspect current state
    print("\n[2/6] Fetching current draft details...")
    draft = make_request(draft_self_url, headers=headers, method="GET")
    if args.draft_id and draft.get("metadata", {}).get("version") not in (
        None,
        "",
        version,
        latest.get("metadata", {}).get("version"),
    ):
        raise ValueError("The selected draft contains a different pending version.")
    today_str = datetime.datetime.now(datetime.timezone.utc).date().isoformat()
    update_payload = build_draft_update(draft, version, today_str)

    # 3. Register files to the draft
    if args.files:
        for fpath in args.files:
            if not os.path.exists(fpath):
                raise FileNotFoundError(f"File not found: {fpath}")
        # Zenodo hands back an unpublished new-version draft left by an
        # earlier run that died later (e.g. a 500 on publish). Files that
        # draft already holds must not be registered again ("already
        # exists"); a stale or half-uploaded one is replaced.
        existing = make_request(draft_files_url, headers=headers, method="GET")
        to_upload, stale = plan_file_uploads(existing.get("entries", []), args.files)
        for key in stale:
            print(f"Removing stale draft file: {key}")
            make_request(
                f"{draft_files_url}/{key}",
                headers=headers,
                method="DELETE",
                json_response=False,
            )
        for fpath in args.files:
            if fpath not in to_upload:
                print(f"Already in the draft, unchanged: {os.path.basename(fpath)}")

        print(f"\n[3/6] Registering {len(to_upload)} files to upload...")
        if to_upload:
            register_payload = [{"key": os.path.basename(f)} for f in to_upload]
            make_request(
                draft_files_url, data=register_payload, headers=headers, method="POST"
            )

        # 4. Upload and commit file contents
        print("\n[4/6] Uploading and committing file contents...")
        for fpath in to_upload:
            if not os.path.exists(fpath):
                raise FileNotFoundError(f"File not found: {fpath}")
            filename = os.path.basename(fpath)
            content_url = f"{draft_files_url}/{filename}/content"
            commit_url = f"{draft_files_url}/{filename}/commit"

            print(f"Uploading file content: {filename} ({content_url})...")
            with open(fpath, "rb") as f:
                file_bytes = f.read()

            # Upload raw file content
            upload_headers = headers.copy()
            upload_headers["Content-Type"] = "application/octet-stream"
            make_request(
                content_url, data=file_bytes, headers=upload_headers, method="PUT"
            )

            # Commit the upload
            print(f"Committing file: {filename}...")
            make_request(commit_url, data=None, headers=headers, method="POST")
    else:
        print("\n[3/6] No files specified to register.")
        print("\n[4/6] No files specified to upload.")

    # 5. Update Metadata (version and publication date)
    print("\n[5/6] Updating metadata...")

    # Read/write the same RDM representation. Rebuilding a draft from the
    # legacy public serializer loses dates, contributor structure and fields.

    print(
        f"Updating draft {draft_id} metadata: version={version}, publication_date={today_str}"
    )
    make_request(draft_self_url, data=update_payload, headers=headers, method="PUT")

    # 6. Publish (Optional)
    if args.publish:
        if not draft_publish_url:
            raise RuntimeError("Publish URL is not available in the draft links.")
        print(f"\n[6/6] Publishing deposition {draft_id}...")
        publish_resp = make_request(
            draft_publish_url, data=None, headers=headers, method="POST"
        )
        final_doi = publish_resp.get("pids", {}).get("doi", {}).get("identifier")
        final_doi = final_doi or publish_resp.get("doi")
        print(f"Successfully published new version! DOI: {final_doi}")
    else:
        print(f"\n[6/6] Deposition {draft_id} updated successfully in draft mode.")


if __name__ == "__main__":
    main()
