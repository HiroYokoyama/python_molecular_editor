# MoleditPy Coverage Report

- **Overall Project Coverage (Full)**: **52.51%**
- **Core Molecular Logic Coverage**: **71.28%**

> [!NOTE]
> **Core Molecular Logic Coverage** excludes UI boilerplate (dialogs, view managers, and interactor styles) to focus on scientific algorithm reliability.

### Core Logic Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\__init__.py             |     18 |      5 |   72.2% |
| moleditpy\src\moleditpy\modules\analysis_window.py      |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\modules\atom_item.py            |    262 |     35 |   86.6% |
| moleditpy\src\moleditpy\modules\bond_item.py            |    340 |     69 |   79.7% |
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    772 |    283 |   63.3% |
| moleditpy\src\moleditpy\modules\constants.py            |     27 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\main_window.py          |    366 |     88 |   76.0% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    396 |    161 |   59.3% |
| moleditpy\src\moleditpy\modules\main_window_compute.py  |    609 |    163 |   73.2% |
| moleditpy\src\moleditpy\modules\main_window_export.py   |    315 |     94 |   70.2% |
| moleditpy\src\moleditpy\modules\main_window_main_init.py |     81 |     38 |   53.1% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    389 |    111 |   71.5% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    173 |     68 |   60.7% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    167 |     32 |   80.8% |
| moleditpy\src\moleditpy\modules\mol_geometry.py         |     92 |     18 |   80.4% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    220 |     39 |   82.3% |
| moleditpy\src\moleditpy\modules\molecule_scene.py       |    461 |    174 |   62.3% |
| moleditpy\src\moleditpy\modules\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\plugin_manager.py       |    282 |     74 |   73.8% |
| **TOTAL** | **5142** | **1477** | **71.28%** |

### Full Application Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\__init__.py             |     18 |      5 |   72.2% |
| moleditpy\src\moleditpy\modules\about_dialog.py         |     65 |     54 |   16.9% |
| moleditpy\src\moleditpy\modules\align_plane_dialog.py   |    161 |     67 |   58.4% |
| moleditpy\src\moleditpy\modules\alignment_dialog.py     |    141 |     51 |   63.8% |
| moleditpy\src\moleditpy\modules\analysis_window.py      |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\modules\angle_dialog.py         |    301 |    165 |   45.2% |
| moleditpy\src\moleditpy\modules\atom_item.py            |    268 |     41 |   84.7% |
| moleditpy\src\moleditpy\modules\bond_item.py            |    352 |     81 |   77.0% |
| moleditpy\src\moleditpy\modules\bond_length_dialog.py   |    272 |    162 |   40.4% |
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    889 |    349 |   60.7% |
| moleditpy\src\moleditpy\modules\color_settings_dialog.py |    250 |    239 |    4.4% |
| moleditpy\src\moleditpy\modules\constants.py            |     27 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\constrained_optimization_dialog.py |    407 |    386 |    5.2% |
| moleditpy\src\moleditpy\modules\custom_interactor_style.py |    478 |    452 |    5.4% |
| moleditpy\src\moleditpy\modules\custom_qt_interactor.py |     43 |     35 |   18.6% |
| moleditpy\src\moleditpy\modules\dialog3_d_picking_mixin.py |    102 |     60 |   41.2% |
| moleditpy\src\moleditpy\modules\dihedral_dialog.py      |    303 |    174 |   42.6% |
| moleditpy\src\moleditpy\modules\main_window.py          |    408 |     96 |   76.5% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    459 |    191 |   58.4% |
| moleditpy\src\moleditpy\modules\main_window_compute.py  |    979 |    459 |   53.1% |
| moleditpy\src\moleditpy\modules\main_window_dialog_manager.py |    243 |    157 |   35.4% |
| moleditpy\src\moleditpy\modules\main_window_edit_3d.py  |    201 |    104 |   48.3% |
| moleditpy\src\moleditpy\modules\main_window_edit_actions.py |    978 |    496 |   49.3% |
| moleditpy\src\moleditpy\modules\main_window_export.py   |    538 |    163 |   69.7% |
| moleditpy\src\moleditpy\modules\main_window_main_init.py |   1374 |    507 |   63.1% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    693 |    374 |   46.0% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    233 |    105 |   54.9% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    167 |     32 |   80.8% |
| moleditpy\src\moleditpy\modules\main_window_ui_manager.py |    326 |    106 |   67.5% |
| moleditpy\src\moleditpy\modules\main_window_view_3d.py  |    926 |    503 |   45.7% |
| moleditpy\src\moleditpy\modules\main_window_view_loaders.py |    211 |    167 |   20.9% |
| moleditpy\src\moleditpy\modules\mirror_dialog.py        |     70 |      7 |   90.0% |
| moleditpy\src\moleditpy\modules\mol_geometry.py         |     92 |     18 |   80.4% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    223 |     42 |   81.2% |
| moleditpy\src\moleditpy\modules\molecule_scene.py       |   1438 |    699 |   51.4% |
| moleditpy\src\moleditpy\modules\move_group_dialog.py    |    384 |    361 |    6.0% |
| moleditpy\src\moleditpy\modules\periodic_table_dialog.py |     33 |     24 |   27.3% |
| moleditpy\src\moleditpy\modules\planarize_dialog.py     |    128 |     37 |   71.1% |
| moleditpy\src\moleditpy\modules\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\plugin_manager.py       |    282 |     74 |   73.8% |
| moleditpy\src\moleditpy\modules\settings_dialog.py      |    819 |    233 |   71.6% |
| moleditpy\src\moleditpy\modules\template_preview_item.py |    101 |     77 |   23.8% |
| moleditpy\src\moleditpy\modules\template_preview_view.py |     41 |     13 |   68.3% |
| moleditpy\src\moleditpy\modules\translation_dialog.py   |    209 |    136 |   34.9% |
| moleditpy\src\moleditpy\modules\user_template_dialog.py |    405 |    167 |   58.8% |
| moleditpy\src\moleditpy\modules\zoomable_view.py        |     72 |     39 |   45.8% |
| **TOTAL** | **16282** | **7733** | **52.51%** |

## Test Suite Status
- **Unit tests**: PASSED
- **Integration tests**: PASSED
- **GUI tests**: PASSED

[View Detailed HTML Report](coverage_html/index.html)