"""Zenodo recovery is idempotent and preserves native RDM draft metadata."""

import copy
import hashlib
import importlib.util
from pathlib import Path
from unittest.mock import MagicMock

import pytest

SPEC = importlib.util.spec_from_file_location(
    "zenodo_recovery", Path(__file__).resolve().parents[2] / "scripts/update_zenodo.py"
)
zen = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(zen)
BASE = "https://zenodo.org/api"


@pytest.fixture
def draft():
    url = BASE + "/records/300/draft"
    return {
        "id": "300",
        "parent": {"id": "family"},
        "metadata": {
            "title": "Synthetic test collection",
            "creators": [{"person_or_org": {"type": "organizational", "name": "Lab"}}],
            "contributors": [
                {"person_or_org": {"name": "Reviewer"}, "role": {"id": "editor"}}
            ],
            "dates": [{"date": "2026-07-11", "type": {"id": "updated"}}],
            "rights": [{"id": "lgpl-3.0-or-later"}],
            "funding": [{"funder": {"name": "Synthetic funder"}}],
            "version": "old",
        },
        "access": {"record": "public", "files": "public", "status": "open"},
        "pids": {
            "doi": {"identifier": "10.5281/zenodo.300", "provider": "datacite"},
            "oai": {"identifier": "oai:zenodo.org:300"},
        },
        "custom_fields": {"code:programmingLanguage": [{"id": "python"}]},
        "links": {
            "self": url,
            "files": url + "/files",
            "publish": url + "/actions/publish",
        },
    }


def test_metadata_preserves_rdm_fields_and_draft_doi(draft):
    before = copy.deepcopy(draft)
    payload = zen.build_draft_update(draft, "next", "2026-10-23")
    expected = dict(before["metadata"], version="next", publication_date="2026-10-23")
    assert payload["metadata"] == expected
    assert payload["pids"] == {"doi": before["pids"]["doi"]}
    assert payload["access"] == {"record": "public", "files": "public"}
    assert payload["custom_fields"] == before["custom_fields"]
    assert set(payload) == {"metadata", "files", "pids", "access", "custom_fields"}
    payload["metadata"]["creators"][0]["person_or_org"]["name"] = "Changed"
    assert draft == before


@pytest.mark.parametrize(
    "right,expected",
    [
        (
            {
                "id": "lgpl-3.0-or-later",
                "title": {"en": "GNU LGPL"},
                "description": {"en": "Vocabulary description"},
                "props": {"scheme": "spdx"},
            },
            {"id": "lgpl-3.0-or-later"},
        ),
        (
            {
                "title": {"en": "Custom license"},
                "description": {"en": "Custom terms"},
                "link": "https://example.org/license",
            },
            {
                "title": {"en": "Custom license"},
                "description": {"en": "Custom terms"},
                "link": "https://example.org/license",
            },
        ),
        ({"id": "gpl-3.0-or-later"}, {"id": "gpl-3.0-or-later"}),
    ],
)
def test_rights_follow_official_id_or_free_text_schema(draft, right, expected):
    draft["metadata"]["rights"] = [copy.deepcopy(right)]
    payload = zen.build_draft_update(draft, "next", "2026-10-23")
    assert payload["metadata"]["rights"] == [expected]
    assert draft["metadata"]["rights"] == [right]


@pytest.mark.parametrize(
    "metadata",
    [
        {},
        {"title": "Missing creators"},
        {"title": "Legacy serializer", "creators": [{"name": "Old format"}]},
    ],
)
def test_incomplete_or_legacy_metadata_is_not_written(metadata):
    with pytest.raises(ValueError):
        zen.build_draft_update({"metadata": metadata}, "next", "2026-10-23")


@pytest.mark.parametrize(
    "path,accept",
    [
        ("/records/1", "application/vnd.inveniordm.v1+json"),
        ("/records/1/draft", "application/vnd.inveniordm.v1+json"),
        ("/records/1/versions/latest", "application/vnd.inveniordm.v1+json"),
        ("/records/1/draft/actions/publish", "application/vnd.inveniordm.v1+json"),
        ("/records/1/draft/files", "application/json"),
        ("/records/1/draft/files/test.zip/content", "application/json"),
    ],
)
def test_native_record_accept_header_and_file_api_header(monkeypatch, path, accept):
    response = MagicMock()
    response.__enter__.return_value.read.return_value = b"{}"
    opener = MagicMock(return_value=response)
    monkeypatch.setattr(zen.urllib.request, "urlopen", opener)
    headers = {"Authorization": "Bearer synthetic-test-token"}
    zen.make_request(BASE + path, headers=headers)
    req = opener.call_args.args[0]
    assert req.get_header("Accept") == accept
    assert headers == {"Authorization": "Bearer synthetic-test-token"}
    assert opener.call_args.kwargs["timeout"] == 60


def invoke(monkeypatch, tmp_path, handler, *args):
    path = tmp_path / "test.zip"
    path.write_bytes(b"synthetic zip")
    monkeypatch.setenv("ZENODO_TOKEN", "synthetic-test-token")
    monkeypatch.setattr(zen, "make_request", handler)
    monkeypatch.setattr(
        "sys.argv", ["upload", "--version-string", "vnext", *args, str(path)]
    )
    zen.main()


@pytest.mark.parametrize("same_asset", [True, False])
def test_already_published_latest_never_creates_a_draft(
    monkeypatch, tmp_path, same_asset
):
    calls = []

    def request(url, **kwargs):
        calls.append((url, kwargs))
        assert url == f"{BASE}/records/{zen.DEPOSITION_ID}/versions/latest"
        checksum = hashlib.md5(
            b"synthetic zip" if same_asset else b"different"
        ).hexdigest()
        return {
            "id": "200",
            "metadata": {"version": "next"},
            "files": {"entries": {"test.zip": {"checksum": "md5:" + checksum}}},
        }

    if same_asset:
        invoke(monkeypatch, tmp_path, request, "--publish")
    else:
        with pytest.raises(ValueError, match="already published"):
            invoke(monkeypatch, tmp_path, request, "--publish")
    assert len(calls) == 1


def test_missing_local_file_does_not_call_remote_api(monkeypatch, tmp_path):
    request = MagicMock()
    monkeypatch.setattr(zen, "make_request", request)
    monkeypatch.setenv("ZENODO_TOKEN", "synthetic-test-token")
    monkeypatch.setattr(
        "sys.argv", ["upload", "--version-string", "next", str(tmp_path / "missing")]
    )
    with pytest.raises(FileNotFoundError):
        zen.main()
    request.assert_not_called()


@pytest.mark.parametrize("same_family", [True, False])
def test_resume_existing_draft_without_creating_a_version(
    monkeypatch, tmp_path, draft, same_family
):
    calls = []
    if not same_family:
        draft["parent"]["id"] = "another-family"

    def request(url, **kwargs):
        method = kwargs.get("method", "GET")
        calls.append((method, url))
        if url.endswith("/versions/latest"):
            return {
                "id": "200",
                "parent": {"id": "family"},
                "metadata": {"version": "old"},
            }
        if method == "GET" and "/versions?" in url:
            return {"hits": {"hits": []}, "links": {}}
        if method == "GET" and url.endswith("/draft"):
            return draft
        if method == "GET" and url.endswith("/files"):
            return {
                "entries": [
                    {
                        "key": "test.zip",
                        "status": "completed",
                        "checksum": "md5:" + hashlib.md5(b"synthetic zip").hexdigest(),
                    }
                ]
            }
        if method == "PUT" and url.endswith("/draft"):
            assert kwargs["data"]["metadata"]["version"] == "next"
            assert kwargs["data"]["pids"]["doi"]["identifier"] == "10.5281/zenodo.300"
            return draft
        if method == "POST" and url.endswith("/actions/publish"):
            return draft
        pytest.fail(f"Unexpected mutation: {method} {url}")

    if same_family:
        invoke(monkeypatch, tmp_path, request, "--draft-id", "300", "--publish")
        assert ("POST", draft["links"]["publish"]) in calls
    else:
        with pytest.raises(ValueError, match="different record family"):
            invoke(monkeypatch, tmp_path, request, "--draft-id", "300", "--publish")
        assert all(method == "GET" for method, _ in calls)
    assert not any(
        method == "POST" and url.endswith("/versions") for method, url in calls
    )


@pytest.mark.parametrize("mapping", [True, False])
def test_existing_file_list_and_mapping_are_both_resumable(tmp_path, mapping):
    path = tmp_path / "test.zip"
    path.write_bytes(b"unchanged")
    entry = {
        "key": path.name,
        "status": "completed",
        "checksum": "md5:" + zen.file_md5(path),
    }
    entries = {path.name: entry} if mapping else [entry]
    assert zen.plan_file_uploads(entries, [str(path)]) == ([], [])


def test_older_published_version_is_found_on_later_page(monkeypatch):
    first = BASE + "/records/200/versions?size=100"
    second = BASE + "/records/200/versions?page=2"
    record = {"id": "100", "metadata": {"version": "older"}}
    responses = {
        first: {"hits": {"hits": []}, "links": {"next": second}},
        second: {"hits": {"hits": [record]}, "links": {}},
        BASE + "/records/100": record,
    }
    calls = []

    def request(url, **kwargs):
        assert kwargs.get("method", "GET") == "GET"
        calls.append(url)
        return responses[url]

    monkeypatch.setattr(zen, "make_request", request)
    assert zen.find_published_version({"id": "200"}, "older", BASE, {}) == record
    assert calls == [first, second, BASE + "/records/100"]


@pytest.mark.parametrize(
    "page",
    [
        {},
        {
            "hits": {"hits": []},
            "links": {"next": BASE + "/records/200/versions?size=100"},
        },
    ],
)
def test_unverifiable_history_fails_before_creating_a_draft(monkeypatch, page):
    monkeypatch.setattr(zen, "make_request", lambda *a, **kw: page)
    with pytest.raises((ValueError, TypeError)):
        zen.find_published_version({"id": "200"}, "next", BASE, {})


def test_new_publication_uses_latest_record_and_uploads_before_publish(
    monkeypatch, tmp_path, draft
):
    calls = []

    def request(url, **kwargs):
        method = kwargs.get("method", "GET")
        calls.append((method, url))
        if url.endswith("/versions/latest"):
            return {"id": "200", "metadata": {"version": "old"}}
        if "/versions?" in url:
            return {"hits": {"hits": []}}
        if method == "POST" and url.endswith("/versions"):
            assert url == BASE + "/records/200/versions"
            return draft
        if method == "GET" and url.endswith("/draft"):
            return draft
        if method == "GET" and url.endswith("/files"):
            return {"entries": []}
        if method == "PUT" and url.endswith("/content"):
            assert kwargs["data"] == b"synthetic zip"
        if method == "POST" and url.endswith("/actions/publish"):
            assert ("POST", draft["links"]["files"] + "/test.zip/commit") in calls
            assert ("PUT", draft["links"]["self"]) in calls
            return draft
        return {}

    invoke(monkeypatch, tmp_path, request, "--publish")
    assert calls[-1] == ("POST", draft["links"]["publish"])
