"""Fresh app imports select the same Qt binding for PyVista and native widgets."""

import os
import subprocess
import sys
import textwrap

import pytest


@pytest.mark.parametrize("previous_binding", [None, "pyqt6"])
def test_app_selects_pyside6_before_pyvistaqt(previous_binding):
    environment = {**os.environ, "QT_QPA_PLATFORM": "offscreen"}
    environment.pop("QT_API", None)
    if previous_binding is not None:
        environment["QT_API"] = previous_binding
    code = textwrap.dedent("""
        import warnings

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            import moleditpy
            import qtpy
            from pyvistaqt import QtInteractor
            from PySide6.QtWidgets import QWidget

        assert qtpy.API_NAME == "PySide6", qtpy.API_NAME
        assert issubclass(QtInteractor, QWidget)
        assert not any(isinstance(item.message, qtpy.PythonQtWarning) for item in caught)
    """)
    result = subprocess.run(
        [sys.executable, "-c", code],
        env=environment,
        capture_output=True,
        text=True,
        errors="replace",
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
