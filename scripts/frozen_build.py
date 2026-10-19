"""Shared frozen-build module restrictions, notices, and reproducibility manifest."""

from importlib import metadata
import json
from pathlib import Path
import pkgutil
import shutil
import subprocess
from typing import Any

import PySide6


def qt_excludes() -> list[str]:
    """Include only the Qt modules used by the application and its renderer."""
    allowed = {
        "QtCore",
        "QtGui",
        "QtWidgets",
        "QtSvg",
        "QtSvgWidgets",
        "QtOpenGL",
        "QtOpenGLWidgets",
        "QtPrintSupport",
        "QtNetwork",
        "QtDBus",
    }
    unused = [
        f"PySide6.{module.name}"
        for module in pkgutil.iter_modules(PySide6.__path__)
        if module.name.startswith("Qt") and module.name not in allowed
    ]
    return ["PyQt6", "PyQt5", "PySide2", "openbabel", *unused]


def license_data(repo_root: Path) -> list[tuple[str, str]]:
    """Preserve application and installed dependency notices in the bundle."""
    output = repo_root / "build" / "distribution-notices"
    output.mkdir(parents=True, exist_ok=True)
    data = []
    for name in (
        "LICENSE",
        "COPYING",
        "LICENSING.md",
        "DISTRIBUTION.md",
        "THIRD_PARTY_NOTICES.md",
    ):
        data.append((str(repo_root / name), "licenses"))
    dependencies: list[dict[str, Any]] = []
    for distribution in sorted(
        metadata.distributions(), key=lambda d: d.metadata["Name"].lower()
    ):
        name = distribution.metadata["Name"]
        dependencies.append(
            {
                "name": name,
                "version": distribution.version,
                "project_urls": distribution.metadata.get_all("Project-URL") or [],
                "license": distribution.metadata.get("License-Expression")
                or distribution.metadata.get("License"),
            }
        )
        for item in distribution.files or []:
            if any(
                word in item.name.lower()
                for word in ("license", "copying", "copyright", "notice")
            ):
                source = Path(distribution.locate_file(item))
                if source.is_file():
                    destination = (
                        output / "third-party" / name / str(item).replace("..", "_")
                    )
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(source, destination)
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=repo_root, text=True
    ).strip()
    manifest = {
        "source_repository": "https://github.com/HiroYokoyama/python_molecular_editor",
        "source_commit": revision,
        "pyside6_version": PySide6.__version__,
        "qt_modules": "Core, Gui, Widgets, Svg and required LGPL support libraries",
        "dependencies": dependencies,
    }
    (output / "build-manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    subprocess.run(
        [
            "git",
            "archive",
            "--format=zip",
            "--output",
            str(output / "moleditpy-source.zip"),
            revision,
        ],
        cwd=repo_root,
        check=True,
    )
    data.append((str(output), "licenses"))
    return data
