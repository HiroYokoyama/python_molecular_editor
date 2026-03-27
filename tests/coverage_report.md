# MoleditPy Coverage Report

- **Overall Project Coverage (Full)**: **63.85%**
- **Core Molecular Logic Coverage**: **80.08%**

> [!NOTE]
> **Core Molecular Logic Coverage** excludes UI boilerplate (dialogs, view managers, and interactor styles) to focus on scientific algorithm reliability.

### Core Logic Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\__init__.py             |     18 |      5 |   72.2% |
| moleditpy\src\moleditpy\modules\analysis_window.py      |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\modules\atom_item.py            |    256 |     33 |   87.1% |
| moleditpy\src\moleditpy\modules\bond_item.py            |    311 |     74 |   76.2% |
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    535 |     86 |   83.9% |
| moleditpy\src\moleditpy\modules\constants.py            |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\main_window.py          |     47 |     18 |   61.7% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    465 |    102 |   78.1% |
| moleditpy\src\moleditpy\modules\main_window_compute.py  |    483 |    112 |   76.8% |
| moleditpy\src\moleditpy\modules\main_window_export.py   |    536 |    161 |   70.0% |
| moleditpy\src\moleditpy\modules\main_window_main_init.py |   1180 |    153 |   87.0% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    350 |     87 |   75.1% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    234 |     79 |   66.2% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    162 |     26 |   84.0% |
| moleditpy\src\moleditpy\modules\mol_geometry.py         |    138 |      7 |   94.9% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    283 |     58 |   79.5% |
| moleditpy\src\moleditpy\modules\molecular_scene_handler.py |    855 |    170 |   80.1% |
| moleditpy\src\moleditpy\modules\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\plugin_manager.py       |    273 |     65 |   76.2% |
| **TOTAL** | **6329** | **1261** | **80.08%** |

### Full Application Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\__init__.py             |     18 |      5 |   72.2% |
| moleditpy\src\moleditpy\modules\about_dialog.py         |     65 |     52 |   20.0% |
| moleditpy\src\moleditpy\modules\align_plane_dialog.py   |    162 |     61 |   62.3% |
| moleditpy\src\moleditpy\modules\alignment_dialog.py     |    142 |     47 |   66.9% |
| moleditpy\src\moleditpy\modules\analysis_window.py      |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\modules\angle_dialog.py         |    298 |    144 |   51.7% |
| moleditpy\src\moleditpy\modules\atom_item.py            |    256 |     33 |   87.1% |
| moleditpy\src\moleditpy\modules\bond_item.py            |    311 |     74 |   76.2% |
| moleditpy\src\moleditpy\modules\bond_length_dialog.py   |    270 |    136 |   49.6% |
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    535 |     86 |   83.9% |
| moleditpy\src\moleditpy\modules\color_settings_dialog.py |    181 |    170 |    6.1% |
| moleditpy\src\moleditpy\modules\constants.py            |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\constrained_optimization_dialog.py |    407 |    307 |   24.6% |
| moleditpy\src\moleditpy\modules\custom_interactor_style.py |    457 |    428 |    6.3% |
| moleditpy\src\moleditpy\modules\custom_qt_interactor.py |     43 |     35 |   18.6% |
| moleditpy\src\moleditpy\modules\dialog_3d_picking_mixin.py |     92 |     50 |   45.7% |
| moleditpy\src\moleditpy\modules\dihedral_dialog.py      |    301 |    139 |   53.8% |
| moleditpy\src\moleditpy\modules\main_window.py          |     47 |     18 |   61.7% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    465 |    102 |   78.1% |
| moleditpy\src\moleditpy\modules\main_window_compute.py  |    483 |    112 |   76.8% |
| moleditpy\src\moleditpy\modules\main_window_dialog_manager.py |    245 |    157 |   35.9% |
| moleditpy\src\moleditpy\modules\main_window_edit_3d.py  |    219 |    121 |   44.7% |
| moleditpy\src\moleditpy\modules\main_window_edit_actions.py |    848 |    366 |   56.8% |
| moleditpy\src\moleditpy\modules\main_window_export.py   |    536 |    161 |   70.0% |
| moleditpy\src\moleditpy\modules\main_window_main_init.py |   1180 |    153 |   87.0% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    350 |     87 |   75.1% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    234 |     79 |   66.2% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    162 |     26 |   84.0% |
| moleditpy\src\moleditpy\modules\main_window_ui_manager.py |    306 |     84 |   72.5% |
| moleditpy\src\moleditpy\modules\main_window_view_3d.py  |    900 |    480 |   46.7% |
| moleditpy\src\moleditpy\modules\main_window_view_loaders.py |    156 |    102 |   34.6% |
| moleditpy\src\moleditpy\modules\mirror_dialog.py        |     70 |      7 |   90.0% |
| moleditpy\src\moleditpy\modules\mol_geometry.py         |    138 |      7 |   94.9% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    283 |     58 |   79.5% |
| moleditpy\src\moleditpy\modules\molecular_scene_handler.py |    855 |    170 |   80.1% |
| moleditpy\src\moleditpy\modules\molecule_scene.py       |    483 |    173 |   64.2% |
| moleditpy\src\moleditpy\modules\move_group_dialog.py    |    377 |    215 |   43.0% |
| moleditpy\src\moleditpy\modules\periodic_table_dialog.py |     33 |      7 |   78.8% |
| moleditpy\src\moleditpy\modules\planarize_dialog.py     |    129 |     37 |   71.3% |
| moleditpy\src\moleditpy\modules\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\plugin_manager.py       |    273 |     65 |   76.2% |
| moleditpy\src\moleditpy\modules\settings_dialog.py      |    712 |    131 |   81.6% |
| moleditpy\src\moleditpy\modules\system_utils.py         |     34 |     18 |   47.1% |
| moleditpy\src\moleditpy\modules\template_preview_item.py |    101 |     77 |   23.8% |
| moleditpy\src\moleditpy\modules\template_preview_view.py |     42 |     14 |   66.7% |
| moleditpy\src\moleditpy\modules\translation_dialog.py   |    197 |     87 |   55.8% |
| moleditpy\src\moleditpy\modules\user_template_dialog.py |    370 |    131 |   64.6% |
| moleditpy\src\moleditpy\modules\zoomable_view.py        |     72 |     39 |   45.8% |
| **TOTAL** | **14041** | **5076** | **63.85%** |

## Test Suite Status
- **Unit tests**: PASSED
- **Integration tests**: PASSED
- **GUI tests**: PASSED

[View Detailed HTML Report](coverage_html/index.html)