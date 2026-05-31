#!/usr/bin/env python
"""
scripts/update_zenodo.py -- Automates updating an existing Zenodo deposition with a new release version
using the new InvenioRDM REST API (/api/records).
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
    parser = argparse.ArgumentParser(description="Update Zenodo Deposition (InvenioRDM API)")
    parser.add_argument("--token", help="Zenodo access token")
    parser.add_argument(
        "--deposition-id",
        default="20463243",
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

    # Normalize version: strip leading v/V if present
    if version.lower().startswith("v"):
        version = version[1:]

    # Determine deposition ID: default to Sandbox ID if --sandbox is used with default production ID
    deposition_id = args.deposition_id
    if args.sandbox and deposition_id == "20463243":
        deposition_id = "506735"

    print(f"Target Version: {version}")
    print(f"Environment: {'Sandbox' if args.sandbox else 'Production'}")
    print(f"Base API: {base_url}")
    print(f"Deposition ID: {deposition_id}")
    print(f"Files to upload: {args.files}")

    headers = {"Authorization": f"Bearer {token}"} if token else {}

    if args.dry_run:
        print("\n--- DRY RUN: Listing Actions ---")
        print(f"[DRY-RUN] Would call POST {base_url}/records/{deposition_id}/versions")
        print(f"[DRY-RUN] Would fetch draft details using GET on the draft link.")
        if args.files:
            print(f"[DRY-RUN] Would register {len(args.files)} files via POST on draft files endpoint.")
            for fpath in args.files:
                print(f"[DRY-RUN] Would upload file: {fpath} via PUT content endpoint.")
                print(f"[DRY-RUN] Would commit file: {fpath} via POST commit endpoint.")
        print(f"[DRY-RUN] Would update metadata with PUT request: version={version}")
        if args.publish:
            print(f"[DRY-RUN] Would publish draft via POST publish endpoint.")
        print("--- DRY RUN COMPLETE ---")
        return

    # 1. Create a new version of the published record
    print(f"\n[1/6] Creating new version draft for record {deposition_id}...")
    new_version_url = f"{base_url}/records/{deposition_id}/versions"
    resp = make_request(new_version_url, data=None, headers=headers, method="POST")

    draft_id = resp.get("id")
    draft_links = resp.get("links", {})
    draft_self_url = draft_links.get("self")
    draft_files_url = draft_links.get("files")
    draft_publish_url = draft_links.get("publish")

    if not draft_id or not draft_self_url:
        raise RuntimeError(f"Could not retrieve draft ID or link from response: {resp}")

    print(f"New Draft ID: {draft_id}")
    print(f"Draft URL: {draft_self_url}")
    print(f"Draft Files URL: {draft_files_url}")

    # 2. Fetch the draft details to inspect current state
    print(f"\n[2/6] Fetching current draft details...")
    draft = make_request(draft_self_url, headers=headers, method="GET")

    # 3. Register files to the draft
    if args.files:
        print(f"\n[3/6] Registering {len(args.files)} files to upload...")
        register_payload = [{"key": os.path.basename(fpath)} for fpath in args.files]
        make_request(draft_files_url, data=register_payload, headers=headers, method="POST")

        # 4. Upload and commit file contents
        print("\n[4/6] Uploading and committing file contents...")
        for fpath in args.files:
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
            make_request(content_url, data=file_bytes, headers=upload_headers, method="PUT")

            # Commit the upload
            print(f"Committing file: {filename}...")
            make_request(commit_url, data=None, headers=headers, method="POST")
    else:
        print("\n[3/6] No files specified to register.")
        print("\n[4/6] No files specified to upload.")

    # 5. Update Metadata (version and publication date)
    print("\n[5/6] Updating metadata...")
    
    # Fetch the parent record to retrieve the original, complete metadata (old info to copy in)
    parent_url = f"{base_url}/records/{deposition_id}"
    print(f"Fetching parent record metadata from {parent_url} to ensure completeness...")
    parent_record = make_request(parent_url, headers=headers, method="GET")
    parent_metadata = parent_record.get("metadata", {})

    # Construct a clean metadata dictionary containing only allowed/editable fields to prevent 500 server errors
    metadata = {}
    
    # Copy essential standard fields from parent record metadata
    for field in ["title", "description", "references"]:
        if field in parent_metadata:
            metadata[field] = parent_metadata[field]

    # Set new version and publisher
    metadata["version"] = version
    metadata["publisher"] = "Zenodo"

    # Set publication date to today's date automatically
    today_str = datetime.date.today().isoformat()
    metadata["publication_date"] = today_str

    # Map and format creators correctly to the InvenioRDM nested structure
    creators = parent_metadata.get("creators", [])
    if isinstance(creators, list):
        new_creators = []
        for c in creators:
            if isinstance(c, dict):
                name = c.get("name")
                if name:
                    c_type = c.get("type", "personal")
                    person_or_org = {
                        "name": name,
                        "type": c_type
                    }
                    if c_type == "personal":
                        if "," in name:
                            parts = name.split(",", 1)
                            family_name = parts[0].strip()
                            given_name = parts[1].strip()
                        elif " " in name:
                            parts = name.rsplit(" ", 1)
                            given_name = parts[0].strip()
                            family_name = parts[1].strip()
                        else:
                            given_name = name
                            family_name = name
                        person_or_org["family_name"] = family_name
                        person_or_org["given_name"] = given_name
                    
                    # Add ORCID and identifiers if present
                    identifiers = []
                    if "orcid" in c and c["orcid"]:
                        identifiers.append({
                            "scheme": "orcid",
                            "identifier": c["orcid"]
                        })
                    if "gnd" in c and c["gnd"]:
                        identifiers.append({
                            "scheme": "gnd",
                            "identifier": c["gnd"]
                        })
                    if identifiers:
                        person_or_org["identifiers"] = identifiers

                    new_creator = {
                        "person_or_org": person_or_org
                    }

                    # Add affiliations if present
                    if "affiliation" in c and c["affiliation"]:
                        new_creator["affiliations"] = [
                            {
                                "name": c["affiliation"]
                            }
                        ]
                    new_creators.append(new_creator)
        metadata["creators"] = new_creators

    # Map dates correctly (each entry needs today's date string and a type dict with an id)
    dates = parent_metadata.get("dates", [])
    if isinstance(dates, list):
        new_dates = []
        for d in dates:
            if isinstance(d, dict):
                d_type = d.get("type")
                
                # Automatically update date in dates list to today's date
                d_date = today_str
                
                # Ensure type is a dictionary with an id
                if isinstance(d_type, str):
                    d_type = {"id": d_type}
                elif isinstance(d_type, dict) and "id" not in d_type:
                    type_id = d_type.get("type") or d_type.get("id") or "other"
                    d_type = {"id": type_id}
                elif not d_type:
                    d_type = {"id": "other"}
                
                new_dates.append({
                    "date": d_date,
                    "type": d_type
                })
        metadata["dates"] = new_dates

    # Map legacy license to InvenioRDM rights list
    lic = parent_metadata.get("license")
    if lic and isinstance(lic, dict):
        lic_id = lic.get("id")
        if lic_id:
            metadata["rights"] = [{"id": lic_id}]

    # Ensure resource_type is in the format expected by the InvenioRDM schema
    resource_type = parent_metadata.get("resource_type", {})
    if isinstance(resource_type, dict):
        rt_type = resource_type.get("type") or resource_type.get("id")
        if rt_type:
            metadata["resource_type"] = {"id": rt_type}

    update_payload = {"metadata": metadata}
    print(f"Updating draft {draft_id} metadata: version={version}, publication_date={today_str}")
    make_request(draft_self_url, data=update_payload, headers=headers, method="PUT")

    # 6. Publish (Optional)
    if args.publish:
        if not draft_publish_url:
            raise RuntimeError("Publish URL is not available in the draft links.")
        print(f"\n[6/6] Publishing deposition {draft_id}...")
        publish_resp = make_request(draft_publish_url, data=None, headers=headers, method="POST")
        final_doi = publish_resp.get("doi")
        print(f"Successfully published new version! DOI: {final_doi}")
    else:
        print(f"\n[6/6] Deposition {draft_id} updated successfully in draft mode.")


if __name__ == "__main__":
    main()
