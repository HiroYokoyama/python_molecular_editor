"""Scene editing must not corrupt interpreter-owned objects in the Qt wheel."""

import os
import subprocess
import sys
import textwrap


def test_native_scene_add_remove_preserves_none_references():
    code = textwrap.dedent("""
        import sys
        from PySide6.QtWidgets import QApplication, QGraphicsScene, QGraphicsRectItem

        app = QApplication([])
        scene = QGraphicsScene()
        item = QGraphicsRectItem()
        for _ in range(10):
            scene.addItem(item)
            scene.removeItem(item)
        before = sys.getrefcount(None)
        for _ in range(200):
            scene.addItem(item)
            scene.removeItem(item)
        after = sys.getrefcount(None)
        assert after >= before, f"Qt scene operations lost {before - after} None references"
    """)
    result = subprocess.run(
        [sys.executable, "-c", code],
        env={**os.environ, "QT_QPA_PLATFORM": "offscreen"},
        capture_output=True,
        text=True,
        errors="replace",
        timeout=30,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
