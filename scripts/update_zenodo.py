#!/usr/bin/env parser
"""
scripts/update_zenodo.py -- Automates updating an existing Zenodo deposition with a new release version.
Supports both Zenodo Sandbox and Production environments, dry-run mode, and metadata updating.
"""
import argparse
import datetime
import json
import os
import re
import urllib.request
import urllib.error


def parse_pyproject_version(filepath="moleditpy/pyproject.toml"):
    """Extract version from pyproject.toml."""
    if not os.path.exists(filepath):
        return None
    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()
    match = re.search(r'(?m)^version\s*=\s*["\']([^"\']+)["\']', content)
    if match:
        return match.group(1)
    return None


def make_request(url, data=None, headers=None, method="GET", json_response=True):
    """Utility function to make HTTP requests using urllib."""
    if headers is None:
        headers = {}

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
        with urllib.request.urlopen(req) as response:
            resp_content = response.read()
            if json_response:
                return json.loads(resp_content.decode("utf-8"))
            return resp_content
    except urllib.error.HTTPError as e:
        err_body = e.read().decode("utf-8", errors="replace")
        raise RuntimeError(
            f"HTTP Error {e.code}: {e.reason}\nURL: {url}\nMethod: {method}\nBody: {err_body}"
        ) from e
    except urllib.error.URLError as e:
        raise RuntimeError(f"URL Error: {e.reason}\nURL: {url}\nMethod: {method}") from e


def main():
    parser = argparse.ArgumentParser(description="Update Zenodo Deposition")
    parser.add_argument("--token", help="Zenodo access token")
    parser.add_argument(
        "--deposition-id",
        default="17268532",
        help="Zenodo deposition ID (concept/latest version ID)",
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
        "https://sandbox.zenodo.org/api"
        if args.sandbox
        else "https://zenodo.org/api"
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

    if version.lower().startswith("v"):
        version = version[1:]

    print(f"Target Version: {version}")
    print(f"Environment: {'Sandbox' if args.sandbox else 'Production'}")
    print(f"Base API: {base_url}")
    print(f"Deposition ID: {args.deposition_id}")
    print(f"Files to upload: {args.files}")

    headers = {"Authorization": f"Bearer {token}"} if token else {}

    if args.dry_run:
        print("\n--- DRY RUN: Listing Actions ---")
        print(f"[DRY-RUN] Would call POST {base_url}/deposit/depositions/{args.deposition_id}/actions/newversion")
        print(f"[DRY-RUN] Would fetch draft details from the returned link.")
        print(f"[DRY-RUN] Would list existing draft files and call DELETE on each of them.")
        for fpath in args.files:
            print(f"[DRY-RUN] Would upload file: {fpath} via PUT request.")
        print(f"[DRY-RUN] Would fetch metadata, update version to {version}, and update publication_date.")
        print(f"[DRY-RUN] Would update metadata with PUT request to the draft deposition endpoint.")
        if args.publish:
            print(f"[DRY-RUN] Would publish the new deposition version via POST action/publish.")
        print("--- DRY RUN COMPLETE ---")
        return

    # 1. Create a new version of the deposition
    print(f"\n[1/5] Creating new version draft for deposition {args.deposition_id}...")
    new_version_url = f"{base_url}/deposit/depositions/{args.deposition_id}/actions/newversion"
    resp = make_request(new_version_url, data={}, headers=headers, method="POST")

    latest_draft_url = resp.get("links", {}).get("latest_draft")
    if not latest_draft_url:
        # Fallback to check if the newversion returned the draft directly
        latest_draft_url = resp.get("links", {}).get("self")
        if "latest_draft" not in resp.get("links", {}):
            raise RuntimeError(
                f"Zenodo response did not contain a latest_draft link. Full response: {resp}"
            )

    # 2. Get the draft deposition details
    print(f"\n[2/5] Fetching details of new version draft: {latest_draft_url}...")
    draft = make_request(latest_draft_url, headers=headers, method="GET")
    draft_id = draft.get("id")
    bucket_url = draft.get("links", {}).get("bucket")

    print(f"New Draft ID: {draft_id}")
    print(f"Bucket URL: {bucket_url}")

    # 3. Clear existing files in the draft deposition
    # In Zenodo, new version drafts sometimes copy files from the previous version.
    # We delete them so we only upload our clean release files.
    existing_files = draft.get("files", [])
    if existing_files:
        print(f"\n[3/5] Cleaning up {len(existing_files)} inherited files from draft...")
        for file_info in existing_files:
            file_id = file_info.get("id")
            # In the old/current API, files can be deleted using their delete link
            delete_url = file_info.get("links", {}).get("self") or f"{latest_draft_url}/files/{file_id}"
            print(f"Deleting file: {file_info.get('filename')} (ID: {file_id})...")
            make_request(delete_url, headers=headers, method="DELETE", json_response=False)
    else:
        print("\n[3/5] No inherited files found in draft to clean up.")

    # 4. Upload files
    if args.files:
        print(f"\n[4/5] Uploading {len(args.files)} release files to the bucket...")
        if not bucket_url:
            raise RuntimeError("Deposition bucket URL is not available for upload.")

        for fpath in args.files:
            if not os.path.exists(fpath):
                raise FileNotFoundError(f"File not found: {fpath}")
            filename = os.path.basename(fpath)
            upload_url = f"{bucket_url}/{filename}"
            print(f"Uploading {filename} to {upload_url}...")
            with open(fpath, "rb") as f:
                file_bytes = f.read()
            # Send raw bytes using PUT request
            make_request(
                upload_url,
                data=file_bytes,
                headers={"Content-Type": "application/octet-stream"},
                method="PUT"
            )
    else:
        print("\n[4/5] No files specified to upload.")

    # 5. Update Metadata
    print("\n[5/5] Updating deposition metadata (version and publication date)...")
    metadata = draft.get("metadata", {}).copy()

    # Update version
    metadata["version"] = version

    # Update publication date to today's date
    today_str = datetime.date.today().isoformat()
    metadata["publication_date"] = today_str

    update_payload = {"metadata": metadata}
    print(f"Updating deposition {draft_id} metadata: version={version}, publication_date={today_str}")
    make_request(latest_draft_url, data=update_payload, headers=headers, method="PUT")

    # 6. Publish (Optional)
    if args.publish:
        print(f"\nPublishing deposition {draft_id}...")
        publish_url = f"{latest_draft_url}/actions/publish"
        publish_resp = make_request(publish_url, data={}, headers=headers, method="POST")
        final_doi = publish_resp.get("doi")
        print(f"Successfully published new version! DOI: {final_doi}")
    else:
        print(f"\nDeposition {draft_id} updated successfully in draft mode.")


if __name__ == "__main__":
    main()
