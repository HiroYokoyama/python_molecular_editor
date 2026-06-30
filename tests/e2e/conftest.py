# -*- coding: utf-8 -*-
"""Conftest for cross-platform integration tests.

CalculationWorker uses only RDKit/Open Babel — no VTK/PyVista mocking needed.
A minimal QApplication fixture is provided for qtbot compatibility.
"""

import os
import sys
import pytest
from PyQt6.QtWidgets import QApplication


@pytest.fixture(scope="session")
def app():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    q_app = QApplication.instance() or QApplication(sys.argv)
    yield q_app
    q_app.quit()
