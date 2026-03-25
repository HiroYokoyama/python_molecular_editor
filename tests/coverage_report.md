# MoleditPy Coverage Report

- **Overall Project Coverage (Full)**: **57.83%**
- **Core Molecular Logic Coverage**: **71.39%**

> [!NOTE]
> **Core Molecular Logic Coverage** excludes UI boilerplate (dialogs, view managers, and interactor styles) to focus on scientific algorithm reliability.

### Core Logic Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\__init__.py             |     18 |      5 |   72.2% |
| moleditpy\src\moleditpy\modules\analysis_window.py      |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\modules\atom_item.py            |    256 |     33 |   87.1% |
| moleditpy\src\moleditpy\modules\bond_item.py            |    308 |     75 |   75.6% |
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    445 |     52 |   88.3% |
| moleditpy\src\moleditpy\modules\constants.py            |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\main_window.py          |    419 |     98 |   76.6% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    450 |    182 |   59.6% |
| moleditpy\src\moleditpy\modules\main_window_compute.py  |    512 |    121 |   76.4% |
| moleditpy\src\moleditpy\modules\main_window_export.py   |    534 |    161 |   69.9% |
| moleditpy\src\moleditpy\modules\main_window_main_init.py |   1291 |    424 |   67.2% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    342 |     79 |   76.9% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    232 |    105 |   54.7% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    160 |     26 |   83.8% |
| moleditpy\src\moleditpy\modules\mol_geometry.py         |     92 |     18 |   80.4% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    272 |     60 |   77.9% |
| moleditpy\src\moleditpy\modules\molecule_scene.py       |   1309 |    505 |   61.4% |
| moleditpy\src\moleditpy\modules\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\plugin_manager.py       |    276 |     68 |   75.4% |
| **TOTAL** | **7119** | **2037** | **71.39%** |

### Full Application Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\__init__.py             |     18 |      5 |   72.2% |
| moleditpy\src\moleditpy\modules\about_dialog.py         |     63 |     52 |   17.5% |
| moleditpy\src\moleditpy\modules\align_plane_dialog.py   |    161 |     67 |   58.4% |
| moleditpy\src\moleditpy\modules\alignment_dialog.py     |    141 |     51 |   63.8% |
| moleditpy\src\moleditpy\modules\analysis_window.py      |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\modules\angle_dialog.py         |    295 |    159 |   46.1% |
| moleditpy\src\moleditpy\modules\atom_item.py            |    256 |     33 |   87.1% |
| moleditpy\src\moleditpy\modules\bond_item.py            |    308 |     75 |   75.6% |
| moleditpy\src\moleditpy\modules\bond_length_dialog.py   |    267 |    157 |   41.2% |
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    445 |     52 |   88.3% |
| moleditpy\src\moleditpy\modules\color_settings_dialog.py |    181 |    170 |    6.1% |
| moleditpy\src\moleditpy\modules\constants.py            |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\constrained_optimization_dialog.py |    405 |    384 |    5.2% |
| moleditpy\src\moleditpy\modules\custom_interactor_style.py |    454 |    428 |    5.7% |
| moleditpy\src\moleditpy\modules\custom_qt_interactor.py |     43 |     35 |   18.6% |
| moleditpy\src\moleditpy\modules\dialog_3d_picking_mixin.py |     92 |     50 |   45.7% |
| moleditpy\src\moleditpy\modules\dihedral_dialog.py      |    298 |    169 |   43.3% |
| moleditpy\src\moleditpy\modules\main_window.py          |    419 |     98 |   76.6% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    450 |    182 |   59.6% |
| moleditpy\src\moleditpy\modules\main_window_compute.py  |    512 |    121 |   76.4% |
| moleditpy\src\moleditpy\modules\main_window_dialog_manager.py |    243 |    157 |   35.4% |
| moleditpy\src\moleditpy\modules\main_window_edit_3d.py  |    217 |    121 |   44.2% |
| moleditpy\src\moleditpy\modules\main_window_edit_actions.py |    855 |    375 |   56.1% |
| moleditpy\src\moleditpy\modules\main_window_export.py   |    534 |    161 |   69.9% |
| moleditpy\src\moleditpy\modules\main_window_main_init.py |   1291 |    424 |   67.2% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    342 |     79 |   76.9% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    232 |    105 |   54.7% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    160 |     26 |   83.8% |
| moleditpy\src\moleditpy\modules\main_window_ui_manager.py |    293 |     79 |   73.0% |
| moleditpy\src\moleditpy\modules\main_window_view_3d.py  |    901 |    481 |   46.6% |
| moleditpy\src\moleditpy\modules\main_window_view_loaders.py |    156 |    102 |   34.6% |
| moleditpy\src\moleditpy\modules\mirror_dialog.py        |     70 |      7 |   90.0% |
| moleditpy\src\moleditpy\modules\mol_geometry.py         |     92 |     18 |   80.4% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    272 |     60 |   77.9% |
| moleditpy\src\moleditpy\modules\molecule_scene.py       |   1309 |    505 |   61.4% |
| moleditpy\src\moleditpy\modules\move_group_dialog.py    |    374 |    351 |    6.1% |
| moleditpy\src\moleditpy\modules\periodic_table_dialog.py |     33 |     24 |   27.3% |
| moleditpy\src\moleditpy\modules\planarize_dialog.py     |    128 |     37 |   71.1% |
| moleditpy\src\moleditpy\modules\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\plugin_manager.py       |    276 |     68 |   75.4% |
| moleditpy\src\moleditpy\modules\settings_dialog.py      |    736 |    157 |   78.7% |
| moleditpy\src\moleditpy\modules\template_preview_item.py |    101 |     77 |   23.8% |
| moleditpy\src\moleditpy\modules\template_preview_view.py |     42 |     14 |   66.7% |
| moleditpy\src\moleditpy\modules\translation_dialog.py   |    196 |    122 |   37.8% |
| moleditpy\src\moleditpy\modules\user_template_dialog.py |    370 |    131 |   64.6% |
| moleditpy\src\moleditpy\modules\zoomable_view.py        |     72 |     39 |   45.8% |
| **TOTAL** | **14306** | **6033** | **57.83%** |

## Test Suite Status
- **Unit tests**: PASSED
- **Integration tests**: PASSED
- **GUI tests**: PASSED

[View Detailed HTML Report](coverage_html/index.html)