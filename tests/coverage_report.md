# MoleditPy Coverage Report

- **Overall Project Coverage (Full)**: **48.82%**
- **Core Molecular Logic Coverage**: **73.21%**

> [!NOTE]
> **Core Molecular Logic Coverage** excludes UI boilerplate (dialogs, view managers, and interactor styles) to focus on scientific algorithm reliability.

### Core Logic Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\__init__.py             |     18 |      5 |   72.2% |
| moleditpy\src\moleditpy\modules\analysis_window.py      |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\modules\atom_item.py            |    261 |     37 |   85.8% |
| moleditpy\src\moleditpy\modules\bond_item.py            |    338 |     64 |   81.1% |
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    478 |    164 |   65.7% |
| moleditpy\src\moleditpy\modules\constants.py            |     27 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\main_window.py          |    362 |     90 |   75.1% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    380 |    145 |   61.8% |
| moleditpy\src\moleditpy\modules\main_window_compute.py  |    529 |     96 |   81.9% |
| moleditpy\src\moleditpy\modules\main_window_export.py   |    315 |     94 |   70.2% |
| moleditpy\src\moleditpy\modules\main_window_main_init.py |     81 |     38 |   53.1% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    389 |    110 |   71.7% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    171 |     66 |   61.4% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    152 |     22 |   85.5% |
| moleditpy\src\moleditpy\modules\mol_geometry.py         |     46 |     16 |   65.2% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    209 |     33 |   84.2% |
| moleditpy\src\moleditpy\modules\molecule_scene.py       |    463 |    176 |   62.0% |
| moleditpy\src\moleditpy\modules\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\plugin_manager.py       |    278 |     70 |   74.8% |
| **TOTAL** | **4669** | **1251** | **73.21%** |

### Full Application Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\__init__.py             |     18 |      5 |   72.2% |
| moleditpy\src\moleditpy\modules\about_dialog.py         |     65 |     54 |   16.9% |
| moleditpy\src\moleditpy\modules\align_plane_dialog.py   |    161 |    141 |   12.4% |
| moleditpy\src\moleditpy\modules\alignment_dialog.py     |    141 |    124 |   12.1% |
| moleditpy\src\moleditpy\modules\analysis_window.py      |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\modules\angle_dialog.py         |    237 |    216 |    8.9% |
| moleditpy\src\moleditpy\modules\atom_item.py            |    267 |     43 |   83.9% |
| moleditpy\src\moleditpy\modules\bond_item.py            |    350 |     76 |   78.3% |
| moleditpy\src\moleditpy\modules\bond_length_dialog.py   |    193 |    175 |    9.3% |
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    536 |    211 |   60.6% |
| moleditpy\src\moleditpy\modules\color_settings_dialog.py |    250 |    239 |    4.4% |
| moleditpy\src\moleditpy\modules\constants.py            |     27 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\constrained_optimization_dialog.py |    405 |    384 |    5.2% |
| moleditpy\src\moleditpy\modules\custom_interactor_style.py |    490 |    475 |    3.1% |
| moleditpy\src\moleditpy\modules\custom_qt_interactor.py |     43 |     35 |   18.6% |
| moleditpy\src\moleditpy\modules\dialog3_d_picking_mixin.py |    102 |     65 |   36.3% |
| moleditpy\src\moleditpy\modules\dihedral_dialog.py      |    223 |    203 |    9.0% |
| moleditpy\src\moleditpy\modules\main_window.py          |    402 |     96 |   76.1% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    443 |    175 |   60.5% |
| moleditpy\src\moleditpy\modules\main_window_compute.py  |    873 |    364 |   58.3% |
| moleditpy\src\moleditpy\modules\main_window_dialog_manager.py |    243 |    157 |   35.4% |
| moleditpy\src\moleditpy\modules\main_window_edit_3d.py  |    234 |    125 |   46.6% |
| moleditpy\src\moleditpy\modules\main_window_edit_actions.py |    973 |    519 |   46.7% |
| moleditpy\src\moleditpy\modules\main_window_export.py   |    536 |    163 |   69.6% |
| moleditpy\src\moleditpy\modules\main_window_main_init.py |   1361 |    504 |   63.0% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    692 |    372 |   46.2% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    231 |    104 |   55.0% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    152 |     22 |   85.5% |
| moleditpy\src\moleditpy\modules\main_window_ui_manager.py |    326 |    117 |   64.1% |
| moleditpy\src\moleditpy\modules\main_window_view_3d.py  |    920 |    669 |   27.3% |
| moleditpy\src\moleditpy\modules\main_window_view_loaders.py |    196 |    153 |   21.9% |
| moleditpy\src\moleditpy\modules\mirror_dialog.py        |     70 |      8 |   88.6% |
| moleditpy\src\moleditpy\modules\mol_geometry.py         |     46 |     16 |   65.2% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    212 |     33 |   84.4% |
| moleditpy\src\moleditpy\modules\molecule_scene.py       |   1437 |    698 |   51.4% |
| moleditpy\src\moleditpy\modules\move_group_dialog.py    |    394 |    371 |    5.8% |
| moleditpy\src\moleditpy\modules\periodic_table_dialog.py |     33 |     24 |   27.3% |
| moleditpy\src\moleditpy\modules\planarize_dialog.py     |    128 |     37 |   71.1% |
| moleditpy\src\moleditpy\modules\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\plugin_manager.py       |    278 |     70 |   74.8% |
| moleditpy\src\moleditpy\modules\settings_dialog.py      |    819 |    233 |   71.6% |
| moleditpy\src\moleditpy\modules\template_preview_item.py |    101 |     77 |   23.8% |
| moleditpy\src\moleditpy\modules\template_preview_view.py |     41 |     13 |   68.3% |
| moleditpy\src\moleditpy\modules\translation_dialog.py   |    223 |    148 |   33.6% |
| moleditpy\src\moleditpy\modules\user_template_dialog.py |    403 |    165 |   59.1% |
| moleditpy\src\moleditpy\modules\zoomable_view.py        |     72 |     39 |   45.8% |
| **TOTAL** | **15519** | **7943** | **48.82%** |

## Test Suite Status
- **Unit tests**: PASSED
- **Integration tests**: PASSED
- **GUI tests**: PASSED

[View Detailed HTML Report](coverage_html/index.html)