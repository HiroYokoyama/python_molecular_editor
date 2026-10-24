"""The catalog must retain tests after component directories are introduced."""

import importlib.util
from pathlib import Path


def test_catalog_finds_nested_tests_in_every_tier(tmp_path: Path) -> None:
    script = Path(__file__).resolve().parents[2] / "utils/generate_assertion_catalog.py"
    spec = importlib.util.spec_from_file_location("assertion_catalog", script)
    assert spec is not None and spec.loader is not None
    catalog = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(catalog)

    tiers = ("unit", "integration", "e2e", "gui", "full_gui")
    for tier in tiers:
        folder = tmp_path / tier / "component"
        folder.mkdir(parents=True)
        (folder / "test_example.py").write_text(
            "class TestExample:\n"
            "    def test_value(self):\n"
            '        """Keep assertions and mock expectations visible."""\n'
            "        assert 1 == 1\n"
            "        callback.assert_called_once()\n",
            encoding="utf-8",
        )
        (folder / "helper.py").write_text(
            "def test_not_collected():\n    assert False\n", encoding="utf-8"
        )

    catalog.generate_catalog(tmp_path)
    report = (tmp_path / "assertion_catalog.md").read_text(encoding="utf-8")
    for tier in tiers:
        assert f"## tests/{tier}/component/test_example.py" in report
    assert report.count("### TestExample.test_value") == len(tiers)
    assert report.count("- assert 1 == 1") == len(tiers)
    assert report.count("- callback.assert_called_once()") == len(tiers)
    assert "test_not_collected" not in report

    catalog.generate_catalog(tmp_path)
    assert (tmp_path / "assertion_catalog.md").read_text(encoding="utf-8") == report
