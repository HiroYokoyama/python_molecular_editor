# MoleditPy Coverage Report

- **Overall Project Coverage (Full)**: **52.79%**
- **Core Molecular Logic Coverage**: **72.00%**

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
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    754 |    265 |   64.9% |
| moleditpy\src\moleditpy\modules\constants.py            |     27 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\main_window.py          |    366 |     88 |   76.0% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    380 |    154 |   59.5% |
| moleditpy\src\moleditpy\modules\main_window_compute.py  |    587 |    154 |   73.8% |
| moleditpy\src\moleditpy\modules\main_window_export.py   |    315 |     94 |   70.2% |
| moleditpy\src\moleditpy\modules\main_window_main_init.py |     81 |     38 |   53.1% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    388 |    110 |   71.6% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    173 |     68 |   60.7% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    152 |     21 |   86.2% |
| moleditpy\src\moleditpy\modules\mol_geometry.py         |     92 |     18 |   80.4% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    213 |     33 |   84.5% |
| moleditpy\src\moleditpy\modules\molecule_scene.py       |    463 |    176 |   62.0% |
| moleditpy\src\moleditpy\modules\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\plugin_manager.py       |    278 |     70 |   74.8% |
| **TOTAL** | **5058** | **1416** | **72.00%** |

### Full Application Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\__init__.py             |     18 |      5 |   72.2% |
| moleditpy\src\moleditpy\modules\about_dialog.py         |     65 |     54 |   16.9% |
| moleditpy\src\moleditpy\modules\align_plane_dialog.py   |    161 |     67 |   58.4% |
| moleditpy\src\moleditpy\modules\alignment_dialog.py     |    141 |     51 |   63.8% |
| moleditpy\src\moleditpy\modules\analysis_window.py      |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\modules\angle_dialog.py         |    296 |    160 |   45.9% |
| moleditpy\src\moleditpy\modules\atom_item.py            |    267 |     40 |   85.0% |
| moleditpy\src\moleditpy\modules\bond_item.py            |    350 |     75 |   78.6% |
| moleditpy\src\moleditpy\modules\bond_length_dialog.py   |    268 |    158 |   41.0% |
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    870 |    330 |   62.1% |
| moleditpy\src\moleditpy\modules\color_settings_dialog.py |    250 |    239 |    4.4% |
| moleditpy\src\moleditpy\modules\constants.py            |     27 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\constrained_optimization_dialog.py |    406 |    385 |    5.2% |
| moleditpy\src\moleditpy\modules\custom_interactor_style.py |    490 |    464 |    5.3% |
| moleditpy\src\moleditpy\modules\custom_qt_interactor.py |     43 |     35 |   18.6% |
| moleditpy\src\moleditpy\modules\dialog3_d_picking_mixin.py |    102 |     60 |   41.2% |
| moleditpy\src\moleditpy\modules\dihedral_dialog.py      |    301 |    172 |   42.9% |
| moleditpy\src\moleditpy\modules\main_window.py          |    406 |     94 |   76.8% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    443 |    184 |   58.5% |
| moleditpy\src\moleditpy\modules\main_window_compute.py  |    945 |    442 |   53.2% |
| moleditpy\src\moleditpy\modules\main_window_dialog_manager.py |    243 |    157 |   35.4% |
| moleditpy\src\moleditpy\modules\main_window_edit_3d.py  |    201 |    104 |   48.3% |
| moleditpy\src\moleditpy\modules\main_window_edit_actions.py |    977 |    495 |   49.3% |
| moleditpy\src\moleditpy\modules\main_window_export.py   |    536 |    163 |   69.6% |
| moleditpy\src\moleditpy\modules\main_window_main_init.py |   1373 |    506 |   63.1% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    691 |    372 |   46.2% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    233 |    105 |   54.9% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    152 |     21 |   86.2% |
| moleditpy\src\moleditpy\modules\main_window_ui_manager.py |    326 |    106 |   67.5% |
| moleditpy\src\moleditpy\modules\main_window_view_3d.py  |    926 |    503 |   45.7% |
| moleditpy\src\moleditpy\modules\main_window_view_loaders.py |    196 |    153 |   21.9% |
| moleditpy\src\moleditpy\modules\mirror_dialog.py        |     70 |      7 |   90.0% |
| moleditpy\src\moleditpy\modules\mol_geometry.py         |     92 |     18 |   80.4% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    216 |     36 |   83.3% |
| moleditpy\src\moleditpy\modules\molecule_scene.py       |   1437 |    698 |   51.4% |
| moleditpy\src\moleditpy\modules\move_group_dialog.py    |    384 |    361 |    6.0% |
| moleditpy\src\moleditpy\modules\periodic_table_dialog.py |     33 |     24 |   27.3% |
| moleditpy\src\moleditpy\modules\planarize_dialog.py     |    128 |     37 |   71.1% |
| moleditpy\src\moleditpy\modules\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\plugin_manager.py       |    278 |     70 |   74.8% |
| moleditpy\src\moleditpy\modules\settings_dialog.py      |    800 |    214 |   73.2% |
| moleditpy\src\moleditpy\modules\template_preview_item.py |    101 |     77 |   23.8% |
| moleditpy\src\moleditpy\modules\template_preview_view.py |     41 |     13 |   68.3% |
| moleditpy\src\moleditpy\modules\translation_dialog.py   |    208 |    135 |   35.1% |
| moleditpy\src\moleditpy\modules\user_template_dialog.py |    403 |    165 |   59.1% |
| moleditpy\src\moleditpy\modules\zoomable_view.py        |     72 |     39 |   45.8% |
| **TOTAL** | **16138** | **7619** | **52.79%** |

## Test Suite Status
- **Unit tests**: PASSED
- **Integration tests**: PASSED
- **GUI tests**: PASSED

[View Detailed HTML Report](coverage_html/index.html)