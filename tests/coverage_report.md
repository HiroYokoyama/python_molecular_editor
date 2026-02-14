# MoleditPy Coverage Report

- **Overall Project Coverage (Full)**: **48.77%**
- **Core Molecular Logic Coverage**: **65.54%**

> [!NOTE]
> **Core Molecular Logic Coverage** excludes UI boilerplate (dialogs, view managers, and interactor styles) to focus on scientific algorithm reliability.

### Core Logic Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\__init__.py             |      8 |      4 |   50.0% |
| moleditpy\src\moleditpy\modules\analysis_window.py      |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\modules\atom_item.py            |    259 |     59 |   77.2% |
| moleditpy\src\moleditpy\modules\bond_item.py            |    341 |    143 |   58.1% |
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    491 |    177 |   64.0% |
| moleditpy\src\moleditpy\modules\constants.py            |     27 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\main_window.py          |    380 |     94 |   75.3% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    408 |    190 |   53.4% |
| moleditpy\src\moleditpy\modules\main_window_compute.py  |    602 |    216 |   64.1% |
| moleditpy\src\moleditpy\modules\main_window_export.py   |    331 |    176 |   46.8% |
| moleditpy\src\moleditpy\modules\main_window_main_init.py |     90 |     46 |   48.9% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    410 |    129 |   68.5% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    179 |     65 |   63.7% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    173 |     31 |   82.1% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    203 |     26 |   87.2% |
| moleditpy\src\moleditpy\modules\molecule_scene.py       |    443 |    173 |   60.9% |
| moleditpy\src\moleditpy\modules\plugin_interface.py     |     58 |     29 |   50.0% |
| moleditpy\src\moleditpy\modules\plugin_manager.py       |    276 |     69 |   75.0% |
| **TOTAL** | **4794** | **1652** | **65.54%** |

### Full Application Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\__init__.py             |     18 |      6 |   66.7% |
| moleditpy\src\moleditpy\modules\about_dialog.py         |     63 |     52 |   17.5% |
| moleditpy\src\moleditpy\modules\align_plane_dialog.py   |    173 |    152 |   12.1% |
| moleditpy\src\moleditpy\modules\alignment_dialog.py     |    157 |    136 |   13.4% |
| moleditpy\src\moleditpy\modules\analysis_window.py      |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\modules\angle_dialog.py         |    275 |    251 |    8.7% |
| moleditpy\src\moleditpy\modules\atom_item.py            |    259 |     59 |   77.2% |
| moleditpy\src\moleditpy\modules\bond_item.py            |    341 |    143 |   58.1% |
| moleditpy\src\moleditpy\modules\bond_length_dialog.py   |    233 |    212 |    9.0% |
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    517 |    195 |   62.3% |
| moleditpy\src\moleditpy\modules\color_settings_dialog.py |    223 |    212 |    4.9% |
| moleditpy\src\moleditpy\modules\constants.py            |     27 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\constrained_optimization_dialog.py |    412 |    390 |    5.3% |
| moleditpy\src\moleditpy\modules\custom_interactor_style.py |    473 |    447 |    5.5% |
| moleditpy\src\moleditpy\modules\custom_qt_interactor.py |     43 |     35 |   18.6% |
| moleditpy\src\moleditpy\modules\dialog3_d_picking_mixin.py |     74 |     53 |   28.4% |
| moleditpy\src\moleditpy\modules\dihedral_dialog.py      |    271 |    249 |    8.1% |
| moleditpy\src\moleditpy\modules\main_window.py          |    420 |     98 |   76.7% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    450 |    199 |   55.8% |
| moleditpy\src\moleditpy\modules\main_window_compute.py  |    821 |    337 |   59.0% |
| moleditpy\src\moleditpy\modules\main_window_dialog_manager.py |    254 |    160 |   37.0% |
| moleditpy\src\moleditpy\modules\main_window_edit_3d.py  |    256 |    125 |   51.2% |
| moleditpy\src\moleditpy\modules\main_window_edit_actions.py |    968 |    509 |   47.4% |
| moleditpy\src\moleditpy\modules\main_window_export.py   |    545 |    253 |   53.6% |
| moleditpy\src\moleditpy\modules\main_window_main_init.py |   1317 |    469 |   64.4% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    679 |    353 |   48.0% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    256 |    117 |   54.3% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    173 |     31 |   82.1% |
| moleditpy\src\moleditpy\modules\main_window_ui_manager.py |    336 |    106 |   68.5% |
| moleditpy\src\moleditpy\modules\main_window_view_3d.py  |    904 |    483 |   46.6% |
| moleditpy\src\moleditpy\modules\main_window_view_loaders.py |    209 |    157 |   24.9% |
| moleditpy\src\moleditpy\modules\mirror_dialog.py        |     65 |      7 |   89.2% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    203 |     26 |   87.2% |
| moleditpy\src\moleditpy\modules\molecule_scene.py       |   1353 |    638 |   52.8% |
| moleditpy\src\moleditpy\modules\move_group_dialog.py    |    381 |    357 |    6.3% |
| moleditpy\src\moleditpy\modules\periodic_table_dialog.py |     33 |     24 |   27.3% |
| moleditpy\src\moleditpy\modules\planarize_dialog.py     |    141 |     38 |   73.0% |
| moleditpy\src\moleditpy\modules\plugin_interface.py     |     58 |     29 |   50.0% |
| moleditpy\src\moleditpy\modules\plugin_manager.py       |    276 |     69 |   75.0% |
| moleditpy\src\moleditpy\modules\settings_dialog.py      |    885 |    297 |   66.4% |
| moleditpy\src\moleditpy\modules\template_preview_item.py |    101 |     77 |   23.8% |
| moleditpy\src\moleditpy\modules\template_preview_view.py |     41 |     13 |   68.3% |
| moleditpy\src\moleditpy\modules\translation_dialog.py   |    225 |    149 |   33.8% |
| moleditpy\src\moleditpy\modules\user_template_dialog.py |    400 |    162 |   59.5% |
| moleditpy\src\moleditpy\modules\zoomable_view.py        |     72 |     39 |   45.8% |
| **TOTAL** | **15496** | **7939** | **48.77%** |

## Test Suite Status
- **Unit tests**: PASSED
- **Integration tests**: PASSED
- **GUI tests**: PASSED

[View Detailed HTML Report](coverage_html/index.html)