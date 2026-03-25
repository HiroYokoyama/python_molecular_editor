# MoleditPy Coverage Report

- **Overall Project Coverage (Full)**: **56.56%**
- **Core Molecular Logic Coverage**: **69.59%**

> [!NOTE]
> **Core Molecular Logic Coverage** excludes UI boilerplate (dialogs, view managers, and interactor styles) to focus on scientific algorithm reliability.

### Core Logic Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\__init__.py             |     18 |      5 |   72.2% |
| moleditpy\src\moleditpy\modules\analysis_window.py      |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\modules\atom_item.py            |    268 |     41 |   84.7% |
| moleditpy\src\moleditpy\modules\bond_item.py            |    352 |     81 |   77.0% |
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    426 |     55 |   87.1% |
| moleditpy\src\moleditpy\modules\constants.py            |     27 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\main_window.py          |    422 |    100 |   76.3% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    456 |    187 |   59.0% |
| moleditpy\src\moleditpy\modules\main_window_compute.py  |    508 |    123 |   75.8% |
| moleditpy\src\moleditpy\modules\main_window_export.py   |    538 |    163 |   69.7% |
| moleditpy\src\moleditpy\modules\main_window_main_init.py |   1330 |    463 |   65.2% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    342 |     87 |   74.6% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    233 |    105 |   54.9% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    167 |     32 |   80.8% |
| moleditpy\src\moleditpy\modules\mol_geometry.py         |     92 |     18 |   80.4% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    217 |     36 |   83.4% |
| moleditpy\src\moleditpy\modules\molecule_scene.py       |   1280 |    573 |   55.2% |
| moleditpy\src\moleditpy\modules\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\plugin_manager.py       |    282 |     74 |   73.8% |
| **TOTAL** | **7130** | **2168** | **69.59%** |

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
| moleditpy\src\moleditpy\modules\atom_item.py            |    268 |     41 |   84.7% |
| moleditpy\src\moleditpy\modules\bond_item.py            |    352 |     81 |   77.0% |
| moleditpy\src\moleditpy\modules\bond_length_dialog.py   |    267 |    157 |   41.2% |
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    426 |     55 |   87.1% |
| moleditpy\src\moleditpy\modules\color_settings_dialog.py |    250 |    239 |    4.4% |
| moleditpy\src\moleditpy\modules\constants.py            |     27 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\constrained_optimization_dialog.py |    405 |    384 |    5.2% |
| moleditpy\src\moleditpy\modules\custom_interactor_style.py |    458 |    432 |    5.7% |
| moleditpy\src\moleditpy\modules\custom_qt_interactor.py |     43 |     35 |   18.6% |
| moleditpy\src\moleditpy\modules\dialog3_d_picking_mixin.py |     92 |     50 |   45.7% |
| moleditpy\src\moleditpy\modules\dihedral_dialog.py      |    298 |    169 |   43.3% |
| moleditpy\src\moleditpy\modules\main_window.py          |    422 |    100 |   76.3% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    456 |    187 |   59.0% |
| moleditpy\src\moleditpy\modules\main_window_compute.py  |    508 |    123 |   75.8% |
| moleditpy\src\moleditpy\modules\main_window_dialog_manager.py |    243 |    157 |   35.4% |
| moleditpy\src\moleditpy\modules\main_window_edit_3d.py  |    217 |    121 |   44.2% |
| moleditpy\src\moleditpy\modules\main_window_edit_actions.py |    854 |    378 |   55.7% |
| moleditpy\src\moleditpy\modules\main_window_export.py   |    538 |    163 |   69.7% |
| moleditpy\src\moleditpy\modules\main_window_main_init.py |   1330 |    463 |   65.2% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    342 |     87 |   74.6% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    233 |    105 |   54.9% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    167 |     32 |   80.8% |
| moleditpy\src\moleditpy\modules\main_window_ui_manager.py |    293 |     80 |   72.7% |
| moleditpy\src\moleditpy\modules\main_window_view_3d.py  |    910 |    489 |   46.3% |
| moleditpy\src\moleditpy\modules\main_window_view_loaders.py |    164 |    111 |   32.3% |
| moleditpy\src\moleditpy\modules\mirror_dialog.py        |     70 |      7 |   90.0% |
| moleditpy\src\moleditpy\modules\mol_geometry.py         |     92 |     18 |   80.4% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    217 |     36 |   83.4% |
| moleditpy\src\moleditpy\modules\molecule_scene.py       |   1280 |    573 |   55.2% |
| moleditpy\src\moleditpy\modules\move_group_dialog.py    |    374 |    351 |    6.1% |
| moleditpy\src\moleditpy\modules\periodic_table_dialog.py |     33 |     24 |   27.3% |
| moleditpy\src\moleditpy\modules\planarize_dialog.py     |    128 |     37 |   71.1% |
| moleditpy\src\moleditpy\modules\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\plugin_manager.py       |    282 |     74 |   73.8% |
| moleditpy\src\moleditpy\modules\settings_dialog.py      |    736 |    157 |   78.7% |
| moleditpy\src\moleditpy\modules\template_preview_item.py |    101 |     77 |   23.8% |
| moleditpy\src\moleditpy\modules\template_preview_view.py |     42 |     14 |   66.7% |
| moleditpy\src\moleditpy\modules\translation_dialog.py   |    196 |    122 |   37.8% |
| moleditpy\src\moleditpy\modules\user_template_dialog.py |    370 |    131 |   64.6% |
| moleditpy\src\moleditpy\modules\zoomable_view.py        |     72 |     39 |   45.8% |
| **TOTAL** | **14406** | **6258** | **56.56%** |

## Test Suite Status
- **Unit tests**: PASSED
- **Integration tests**: PASSED
- **GUI tests**: PASSED

[View Detailed HTML Report](coverage_html/index.html)