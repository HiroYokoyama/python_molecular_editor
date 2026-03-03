# MoleditPy Coverage Report

- **Overall Project Coverage (Full)**: **50.92%**
- **Core Molecular Logic Coverage**: **72.55%**

> [!NOTE]
> **Core Molecular Logic Coverage** excludes UI boilerplate (dialogs, view managers, and interactor styles) to focus on scientific algorithm reliability.

### Core Logic Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\__init__.py             |     18 |      5 |   72.2% |
| moleditpy\src\moleditpy\modules\analysis_window.py      |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\modules\atom_item.py            |    261 |     34 |   87.0% |
| moleditpy\src\moleditpy\modules\bond_item.py            |    338 |     63 |   81.4% |
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    739 |    237 |   67.9% |
| moleditpy\src\moleditpy\modules\constants.py            |     27 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\main_window.py          |    364 |     87 |   76.1% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    380 |    145 |   61.8% |
| moleditpy\src\moleditpy\modules\main_window_compute.py  |    587 |    155 |   73.6% |
| moleditpy\src\moleditpy\modules\main_window_export.py   |    315 |     94 |   70.2% |
| moleditpy\src\moleditpy\modules\main_window_main_init.py |     81 |     38 |   53.1% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    388 |    110 |   71.6% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    173 |     68 |   60.7% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    152 |     21 |   86.2% |
| moleditpy\src\moleditpy\modules\mol_geometry.py         |     79 |     17 |   78.5% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    212 |     35 |   83.5% |
| moleditpy\src\moleditpy\modules\molecule_scene.py       |    463 |    176 |   62.0% |
| moleditpy\src\moleditpy\modules\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\plugin_manager.py       |    278 |     70 |   74.8% |
| **TOTAL** | **5027** | **1380** | **72.55%** |

### Full Application Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\__init__.py             |     18 |      5 |   72.2% |
| moleditpy\src\moleditpy\modules\about_dialog.py         |     65 |     54 |   16.9% |
| moleditpy\src\moleditpy\modules\align_plane_dialog.py   |    161 |    141 |   12.4% |
| moleditpy\src\moleditpy\modules\alignment_dialog.py     |    141 |    124 |   12.1% |
| moleditpy\src\moleditpy\modules\analysis_window.py      |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\modules\angle_dialog.py         |    298 |    190 |   36.2% |
| moleditpy\src\moleditpy\modules\atom_item.py            |    267 |     40 |   85.0% |
| moleditpy\src\moleditpy\modules\bond_item.py            |    350 |     75 |   78.6% |
| moleditpy\src\moleditpy\modules\bond_length_dialog.py   |    274 |    250 |    8.8% |
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    844 |    308 |   63.5% |
| moleditpy\src\moleditpy\modules\color_settings_dialog.py |    250 |    239 |    4.4% |
| moleditpy\src\moleditpy\modules\constants.py            |     27 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\constrained_optimization_dialog.py |    405 |    384 |    5.2% |
| moleditpy\src\moleditpy\modules\custom_interactor_style.py |    490 |    464 |    5.3% |
| moleditpy\src\moleditpy\modules\custom_qt_interactor.py |     43 |     35 |   18.6% |
| moleditpy\src\moleditpy\modules\dialog3_d_picking_mixin.py |    102 |     65 |   36.3% |
| moleditpy\src\moleditpy\modules\dihedral_dialog.py      |    303 |    195 |   35.6% |
| moleditpy\src\moleditpy\modules\main_window.py          |    404 |     93 |   77.0% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    443 |    175 |   60.5% |
| moleditpy\src\moleditpy\modules\main_window_compute.py  |    933 |    432 |   53.7% |
| moleditpy\src\moleditpy\modules\main_window_dialog_manager.py |    243 |    157 |   35.4% |
| moleditpy\src\moleditpy\modules\main_window_edit_3d.py  |    234 |    125 |   46.6% |
| moleditpy\src\moleditpy\modules\main_window_edit_actions.py |    977 |    495 |   49.3% |
| moleditpy\src\moleditpy\modules\main_window_export.py   |    536 |    163 |   69.6% |
| moleditpy\src\moleditpy\modules\main_window_main_init.py |   1363 |    503 |   63.1% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    691 |    372 |   46.2% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    233 |    105 |   54.9% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    152 |     21 |   86.2% |
| moleditpy\src\moleditpy\modules\main_window_ui_manager.py |    326 |    106 |   67.5% |
| moleditpy\src\moleditpy\modules\main_window_view_3d.py  |    926 |    503 |   45.7% |
| moleditpy\src\moleditpy\modules\main_window_view_loaders.py |    196 |    153 |   21.9% |
| moleditpy\src\moleditpy\modules\mirror_dialog.py        |     70 |      7 |   90.0% |
| moleditpy\src\moleditpy\modules\mol_geometry.py         |     79 |     17 |   78.5% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    215 |     35 |   83.7% |
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
| **TOTAL** | **16160** | **7931** | **50.92%** |

## Test Suite Status
- **Unit tests**: PASSED
- **Integration tests**: PASSED
- **GUI tests**: PASSED

[View Detailed HTML Report](coverage_html/index.html)