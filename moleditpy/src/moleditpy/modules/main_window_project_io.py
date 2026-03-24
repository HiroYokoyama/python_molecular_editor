#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

"""
main_window_project_io.py
Functional class separated from main_window.py
"""

import copy
import json
import os
import pickle
import traceback

# PyQt6 Modules
from PyQt6.QtCore import QTimer
from PyQt6.QtWidgets import QFileDialog, QMessageBox

try:
    from PyQt6 import sip as _sip  # type: ignore
    _sip_isdeleted = getattr(_sip, "isdeleted", None)
except Exception:
    _sip = None
    _sip_isdeleted = None


# --- Classes ---
class MainWindowProjectIo(object):
    """Functional class separated from main_window.py."""

    def save_project(self):
        """Save (Ctrl+S) - Defaults to PMEPRJ format."""
        if not self.data.atoms and not self.current_mol:
            self.statusBar().showMessage("Error: Nothing to save.")
            return
        # Non-native formats (MOL/SDF/XYZ) trigger "Save As"
        native_exts = [".pmeprj", ".pmeraw"]
        if self.current_file_path and any(
            self.current_file_path.lower().endswith(ext) for ext in native_exts
        ):
            try:
                if self.current_file_path.lower().endswith(".pmeraw"):
                    save_data = self.get_current_state()
                    with open(self.current_file_path, "wb") as f:
                        pickle.dump(save_data, f)
                else:
                    # Save in PMEPRJ format
                    json_data = self.create_json_data()
                    with open(self.current_file_path, "w", encoding="utf-8") as f:
                        json.dump(json_data, f, indent=2, ensure_ascii=False)

                # Reset state on success
                self.has_unsaved_changes = False
                self.update_window_title()

                self.statusBar().showMessage(
                    f"Project saved to {self.current_file_path}"
                )

            except (OSError, IOError) as e:  # pragma: no cover
                self.statusBar().showMessage(f"File I/O error: {e}")
            except (
                pickle.PicklingError,
                TypeError,
                ValueError,
            ) as e:  # pragma: no cover
                self.statusBar().showMessage(f"Data serialization error: {e}")
            except Exception as e:
                self.statusBar().showMessage(f"Error saving project file: {e}")
        else:
            # Force "Save As" for non-native formats
            self.save_project_as()

    def save_project_as(self):
        """Save As (Ctrl+Shift+S) - Defaults to PMEPRJ format."""
        if not self.data.atoms and not self.current_mol:
            self.statusBar().showMessage("Error: Nothing to save.")
            return

        try:
            # Determine a sensible default filename based on current file (strip extension)
            default_name = "untitled"
            try:
                if self.current_file_path:
                    base = os.path.basename(self.current_file_path)
                    default_name = os.path.splitext(base)[0]
            except Exception:
                default_name = "untitled"

            # Prefer the directory of the currently opened file as default
            default_path = default_name
            try:
                if self.current_file_path:
                    default_path = os.path.join(
                        os.path.dirname(self.current_file_path), default_name
                    )
            except Exception:
                default_path = default_name

            file_path, _ = QFileDialog.getSaveFileName(  # pragma: no cover
                self,
                "Save Project As",
                default_path,
                "PME Project Files (*.pmeprj);;All Files (*)",
            )
            if not file_path:  # pragma: no cover
                return

            if not file_path.lower().endswith(".pmeprj"):
                file_path += ".pmeprj"

            json_data = self.create_json_data()
            with open(file_path, "w", encoding="utf-8") as f:
                json.dump(json_data, f, indent=2, ensure_ascii=False)

            # Reset state on success
            self.has_unsaved_changes = False
            # Replace current file with the newly saved file so subsequent saves go to this path
            self.current_file_path = file_path
            self.update_window_title()
            # Mark this state as the last saved state for undo tracking
            try:
                self._saved_state = copy.deepcopy(self.get_current_state())
            except Exception:  # pragma: no cover
                traceback.print_exc()
            self.statusBar().showMessage(f"Project saved to {file_path}")

        except (OSError, IOError) as e:  # pragma: no cover
            self.statusBar().showMessage(f"File I/O error: {e}")
        except pickle.PicklingError as e:  # pragma: no cover
            self.statusBar().showMessage(f"Data serialization error: {e}")
        except Exception as e:
            self.statusBar().showMessage(f"Error saving project file: {e}")

    def save_raw_data(self):
        if not self.data.atoms and not self.current_mol:
            self.statusBar().showMessage("Error: Nothing to save.")
            return

        try:
            save_data = self.get_current_state()
            # default filename based on current file
            default_name = "untitled"
            try:
                if self.current_file_path:
                    base = os.path.basename(self.current_file_path)
                    default_name = os.path.splitext(base)[0]
            except Exception:
                default_name = "untitled"

            # prefer same directory as current file when available
            default_path = default_name
            try:
                if self.current_file_path:
                    default_path = os.path.join(
                        os.path.dirname(self.current_file_path), default_name
                    )
            except Exception:
                default_path = default_name

            file_path, _ = QFileDialog.getSaveFileName(
                self,
                "Save Project File",
                default_path,
                "Project Files (*.pmeraw);;All Files (*)",
            )  # pragma: no cover
            if not file_path:  # pragma: no cover
                return

            if not file_path.lower().endswith(".pmeraw"):
                file_path += ".pmeraw"

            with open(file_path, "wb") as f:
                pickle.dump(save_data, f)

            # Reset state on success
            self.has_unsaved_changes = False
            # Update current file to the newly saved raw file
            self.current_file_path = file_path
            self.update_window_title()
            try:
                self._saved_state = copy.deepcopy(self.get_current_state())
            except Exception:  # pragma: no cover
                traceback.print_exc()
            self.statusBar().showMessage(f"Project saved to {file_path}")

        except (OSError, IOError) as e:  # pragma: no cover
            self.statusBar().showMessage(f"File I/O error: {e}")
        except pickle.PicklingError as e:  # pragma: no cover
            self.statusBar().showMessage(f"Data serialization error: {e}")
        except Exception as e:
            self.statusBar().showMessage(f"Error saving project file: {e}")

    def load_raw_data(self, file_path=None):
        if not file_path:  # pragma: no cover
            file_path, _ = QFileDialog.getOpenFileName(
                self, "Open Project File", "", "Project Files (*.pmeraw);;All Files (*)"
            )
            if not file_path:
                return

        if not self.clear_all():
            return

        try:
            with open(file_path, "rb") as f:
                loaded_data = pickle.load(f)
            self.restore_ui_for_editing()
            self.set_state_from_data(loaded_data)

            # Reset state on load
            self.reset_undo_stack()
            self.has_unsaved_changes = False
            self.current_file_path = file_path
            self.update_window_title()
            try:
                self._saved_state = copy.deepcopy(self.get_current_state())
            except Exception:  # pragma: no cover
                traceback.print_exc()
            self.statusBar().showMessage(f"Project loaded from {file_path}")

            QTimer.singleShot(0, self.fit_to_view)

        except FileNotFoundError:  # pragma: no cover
            self.statusBar().showMessage(f"File not found: {file_path}")
        except (OSError, IOError) as e:  # pragma: no cover
            self.statusBar().showMessage(f"File I/O error: {e}")
        except pickle.UnpicklingError as e:  # pragma: no cover
            self.statusBar().showMessage(f"Invalid project file format: {e}")
        except Exception as e:
            self.statusBar().showMessage(f"Error loading project file: {e}")

    def save_as_json(self):
        """Save as PMEProject (JSON) with 3D info."""
        if not self.data.atoms and not self.current_mol:
            self.statusBar().showMessage("Error: Nothing to save.")
            return

        try:
            # default filename based on current file
            default_name = "untitled"
            try:
                if self.current_file_path:
                    base = os.path.basename(self.current_file_path)
                    default_name = os.path.splitext(base)[0]
            except Exception:
                default_name = "untitled"

            # prefer same directory as current file when available
            default_path = default_name
            try:
                if self.current_file_path:
                    default_path = os.path.join(
                        os.path.dirname(self.current_file_path), default_name
                    )
            except Exception:
                default_path = default_name

            file_path, _ = QFileDialog.getSaveFileName(  # pragma: no cover
                self,
                "Save as PME Project",
                default_path,
                "PME Project Files (*.pmeprj);;All Files (*)",
            )
            if not file_path:  # pragma: no cover
                return

            if not file_path.lower().endswith(".pmeprj"):
                file_path += ".pmeprj"

            json_data = self.create_json_data()

            with open(file_path, "w", encoding="utf-8") as f:
                json.dump(json_data, f, indent=2, ensure_ascii=False)

            # Reset state on success
            self.has_unsaved_changes = False
            # Replace current file with the newly saved PME Project
            self.current_file_path = file_path
            self.update_window_title()

            self.statusBar().showMessage(f"PME Project saved to {file_path}")

        except (OSError, IOError) as e:  # pragma: no cover
            self.statusBar().showMessage(f"File I/O error: {e}")
        except (TypeError, ValueError) as e:  # pragma: no cover
            self.statusBar().showMessage(f"JSON serialization error: {e}")
        except Exception as e:
            self.statusBar().showMessage(f"Error saving PME Project file: {e}")

    def load_json_data(self, file_path=None):
        """Load PME Project file."""
        if not file_path:  # pragma: no cover
            file_path, _ = QFileDialog.getOpenFileName(
                self,
                "Open PME Project File",
                "",
                "PME Project Files (*.pmeprj);;All Files (*)",
            )
            if not file_path:
                return

        if not self.clear_all():
            return

        try:
            with open(file_path, "r", encoding="utf-8") as f:
                json_data = json.load(f)

            # Format validation
            if json_data.get("format") != "PME Project":  # pragma: no cover
                QMessageBox.warning(
                    self,
                    "Invalid Format",
                    "This file is not a valid PME Project format.",
                )
                return

            # Version check
            file_version = json_data.get("version", "1.0")
            if file_version != "1.0":  # pragma: no cover
                QMessageBox.information(
                    self,
                    "Version Notice",
                    f"This file was created with PME Project version {file_version}.\n"
                    "Loading will be attempted but some features may not work correctly.",
                )

            self.restore_ui_for_editing()
            self.load_from_json_data(json_data)
            # Reset state on load
            self.reset_undo_stack()
            self.has_unsaved_changes = False
            self.current_file_path = file_path
            self.update_window_title()

            self.statusBar().showMessage(f"PME Project loaded from {file_path}")

            QTimer.singleShot(0, self.fit_to_view)

        except FileNotFoundError:  # pragma: no cover
            self.statusBar().showMessage(f"File not found: {file_path}")
        except json.JSONDecodeError as e:  # pragma: no cover
            self.statusBar().showMessage(f"Invalid JSON format: {e}")
        except (OSError, IOError) as e:  # pragma: no cover
            self.statusBar().showMessage(f"File I/O error: {e}")
        except Exception as e:
            self.statusBar().showMessage(f"Error loading PME Project file: {e}")

    def open_project_file(self, file_path=None):
        """Open project file (.pmeprj or .pmeraw)."""
        if not file_path:  # pragma: no cover
            file_path, _ = QFileDialog.getOpenFileName(
                self,
                "Open Project File",
                "",
                "PME Project Files (*.pmeprj);;PME Raw Files (*.pmeraw);;All Files (*)",
            )
            if not file_path:
                return

        # Dispatch based on extension
        if file_path.lower().endswith(".pmeprj"):
            self.load_json_data(file_path)
        elif file_path.lower().endswith(".pmeraw"):
            self.load_raw_data(file_path)
        else:
            # Try JSON if extension unknown
            try:
                self.load_json_data(file_path)
            except Exception:
                try:
                    self.load_raw_data(file_path)
                except Exception:
                    self.statusBar().showMessage(
                        "Error: Unable to determine file format."
                    )
