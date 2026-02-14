# MoleditPy Coverage Report

- **Overall Project Coverage (Full)**: **50.06%**
- **Core Molecular Logic Coverage**: **71.03%**

> [!NOTE]
> **Core Molecular Logic Coverage** excludes UI boilerplate (dialogs, view managers, and interactor styles) to focus on scientific algorithm reliability.

### Core Logic Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\__init__.py             |      8 |      4 |   50.0% |
| moleditpy\src\moleditpy\modules\analysis_window.py      |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\modules\atom_item.py            |    259 |     38 |   85.3% |
| moleditpy\src\moleditpy\modules\bond_item.py            |    340 |     71 |   79.1% |
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    491 |    176 |   64.2% |
| moleditpy\src\moleditpy\modules\constants.py            |     27 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\main_window.py          |    366 |     90 |   75.4% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    396 |    161 |   59.3% |
| moleditpy\src\moleditpy\modules\main_window_compute.py  |    581 |    172 |   70.4% |
| moleditpy\src\moleditpy\modules\main_window_export.py   |    319 |    105 |   67.1% |
| moleditpy\src\moleditpy\modules\main_window_main_init.py |     83 |     41 |   50.6% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    397 |    122 |   69.3% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    168 |     58 |   65.5% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    161 |     25 |   84.5% |
| moleditpy\src\moleditpy\modules\mol_geometry.py         |     46 |     16 |   65.2% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    203 |     26 |   87.2% |
| moleditpy\src\moleditpy\modules\molecule_scene.py       |    443 |    173 |   60.9% |
| moleditpy\src\moleditpy\modules\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\plugin_manager.py       |    276 |     69 |   75.0% |
| **TOTAL** | **4736** | **1372** | **71.03%** |

### Full Application Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\__init__.py             |     18 |      6 |   66.7% |
| moleditpy\src\moleditpy\modules\about_dialog.py         |     63 |     52 |   17.5% |
| moleditpy\src\moleditpy\modules\align_plane_dialog.py   |    160 |    140 |   12.5% |
| moleditpy\src\moleditpy\modules\alignment_dialog.py     |    142 |    124 |   12.7% |
| moleditpy\src\moleditpy\modules\analysis_window.py      |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\modules\angle_dialog.py         |    233 |    212 |    9.0% |
| moleditpy\src\moleditpy\modules\atom_item.py            |    259 |     38 |   85.3% |
| moleditpy\src\moleditpy\modules\bond_item.py            |    340 |     71 |   79.1% |
| moleditpy\src\moleditpy\modules\bond_length_dialog.py   |    190 |    172 |    9.5% |
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    517 |    194 |   62.5% |
| moleditpy\src\moleditpy\modules\color_settings_dialog.py |    225 |    214 |    4.9% |
| moleditpy\src\moleditpy\modules\constants.py            |     27 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\constrained_optimization_dialog.py |    404 |    383 |    5.2% |
| moleditpy\src\moleditpy\modules\custom_interactor_style.py |    473 |    458 |    3.2% |
| moleditpy\src\moleditpy\modules\custom_qt_interactor.py |     43 |     35 |   18.6% |
| moleditpy\src\moleditpy\modules\dialog3_d_picking_mixin.py |     97 |     60 |   38.1% |
| moleditpy\src\moleditpy\modules\dihedral_dialog.py      |    220 |    200 |    9.1% |
| moleditpy\src\moleditpy\modules\main_window.py          |    406 |     94 |   76.8% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    438 |    170 |   61.2% |
| moleditpy\src\moleditpy\modules\main_window_compute.py  |    800 |    291 |   63.6% |
| moleditpy\src\moleditpy\modules\main_window_dialog_manager.py |    242 |    153 |   36.8% |
| moleditpy\src\moleditpy\modules\main_window_edit_3d.py  |    228 |    119 |   47.8% |
| moleditpy\src\moleditpy\modules\main_window_edit_actions.py |    949 |    493 |   48.1% |
| moleditpy\src\moleditpy\modules\main_window_export.py   |    533 |    172 |   67.7% |
| moleditpy\src\moleditpy\modules\main_window_main_init.py |   1310 |    464 |   64.6% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    666 |    346 |   48.0% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    245 |    110 |   55.1% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    161 |     25 |   84.5% |
| moleditpy\src\moleditpy\modules\main_window_ui_manager.py |    324 |    110 |   66.0% |
| moleditpy\src\moleditpy\modules\main_window_view_3d.py  |    892 |    639 |   28.4% |
| moleditpy\src\moleditpy\modules\main_window_view_loaders.py |    197 |    150 |   23.9% |
| moleditpy\src\moleditpy\modules\mirror_dialog.py        |     65 |      8 |   87.7% |
| moleditpy\src\moleditpy\modules\mol_geometry.py         |     46 |     16 |   65.2% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    203 |     26 |   87.2% |
| moleditpy\src\moleditpy\modules\molecule_scene.py       |   1353 |    638 |   52.8% |
| moleditpy\src\moleditpy\modules\move_group_dialog.py    |    381 |    358 |    6.0% |
| moleditpy\src\moleditpy\modules\periodic_table_dialog.py |     33 |     24 |   27.3% |
| moleditpy\src\moleditpy\modules\planarize_dialog.py     |    128 |     37 |   71.1% |
| moleditpy\src\moleditpy\modules\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\plugin_manager.py       |    276 |     69 |   75.0% |
| moleditpy\src\moleditpy\modules\settings_dialog.py      |    787 |    200 |   74.6% |
| moleditpy\src\moleditpy\modules\template_preview_item.py |    101 |     77 |   23.8% |
| moleditpy\src\moleditpy\modules\template_preview_view.py |     41 |     13 |   68.3% |
| moleditpy\src\moleditpy\modules\translation_dialog.py   |    219 |    143 |   34.7% |
| moleditpy\src\moleditpy\modules\user_template_dialog.py |    400 |    162 |   59.5% |
| moleditpy\src\moleditpy\modules\zoomable_view.py        |     72 |     39 |   45.8% |
| **TOTAL** | **15079** | **7530** | **50.06%** |

## Test Suite Status
- **Unit tests**: PASSED
- **Integration tests**: PASSED
- **GUI tests**: PASSED

[View Detailed HTML Report](coverage_html/index.html)