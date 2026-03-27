# MoleditPy Coverage Report

- **Overall Project Coverage (Full)**: **60.51%**
- **Core Molecular Logic Coverage**: **78.37%**

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
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    470 |     57 |   87.9% |
| moleditpy\src\moleditpy\modules\constants.py            |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\main_window.py          |     47 |     18 |   61.7% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    468 |    108 |   76.9% |
| moleditpy\src\moleditpy\modules\main_window_compute.py  |    482 |    110 |   77.2% |
| moleditpy\src\moleditpy\modules\main_window_export.py   |    536 |    161 |   70.0% |
| moleditpy\src\moleditpy\modules\main_window_main_init.py |   1170 |    148 |   87.4% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    343 |     79 |   77.0% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    234 |    105 |   55.1% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    162 |     26 |   84.0% |
| moleditpy\src\moleditpy\modules\mol_geometry.py         |    138 |     19 |   86.2% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    273 |     55 |   79.9% |
| moleditpy\src\moleditpy\modules\molecular_scene_handler.py |    846 |    257 |   69.6% |
| moleditpy\src\moleditpy\modules\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\plugin_manager.py       |    276 |     68 |   75.4% |
| **TOTAL** | **6233** | **1348** | **78.37%** |

### Full Application Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\__init__.py             |     18 |      5 |   72.2% |
| moleditpy\src\moleditpy\modules\about_dialog.py         |     65 |     52 |   20.0% |
| moleditpy\src\moleditpy\modules\align_plane_dialog.py   |    161 |     67 |   58.4% |
| moleditpy\src\moleditpy\modules\alignment_dialog.py     |    141 |     51 |   63.8% |
| moleditpy\src\moleditpy\modules\analysis_window.py      |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\modules\angle_dialog.py         |    297 |    159 |   46.5% |
| moleditpy\src\moleditpy\modules\atom_item.py            |    256 |     33 |   87.1% |
| moleditpy\src\moleditpy\modules\bond_item.py            |    311 |     74 |   76.2% |
| moleditpy\src\moleditpy\modules\bond_length_dialog.py   |    269 |    157 |   41.6% |
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    470 |     57 |   87.9% |
| moleditpy\src\moleditpy\modules\color_settings_dialog.py |    181 |    170 |    6.1% |
| moleditpy\src\moleditpy\modules\constants.py            |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\constrained_optimization_dialog.py |    407 |    384 |    5.7% |
| moleditpy\src\moleditpy\modules\custom_interactor_style.py |    456 |    428 |    6.1% |
| moleditpy\src\moleditpy\modules\custom_qt_interactor.py |     43 |     35 |   18.6% |
| moleditpy\src\moleditpy\modules\dialog_3d_picking_mixin.py |     92 |     50 |   45.7% |
| moleditpy\src\moleditpy\modules\dihedral_dialog.py      |    300 |    169 |   43.7% |
| moleditpy\src\moleditpy\modules\main_window.py          |     47 |     18 |   61.7% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    468 |    108 |   76.9% |
| moleditpy\src\moleditpy\modules\main_window_compute.py  |    482 |    110 |   77.2% |
| moleditpy\src\moleditpy\modules\main_window_dialog_manager.py |    245 |    157 |   35.9% |
| moleditpy\src\moleditpy\modules\main_window_edit_3d.py  |    219 |    121 |   44.7% |
| moleditpy\src\moleditpy\modules\main_window_edit_actions.py |    847 |    366 |   56.8% |
| moleditpy\src\moleditpy\modules\main_window_export.py   |    536 |    161 |   70.0% |
| moleditpy\src\moleditpy\modules\main_window_main_init.py |   1170 |    148 |   87.4% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    343 |     79 |   77.0% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    234 |    105 |   55.1% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    162 |     26 |   84.0% |
| moleditpy\src\moleditpy\modules\main_window_ui_manager.py |    304 |     82 |   73.0% |
| moleditpy\src\moleditpy\modules\main_window_view_3d.py  |    900 |    480 |   46.7% |
| moleditpy\src\moleditpy\modules\main_window_view_loaders.py |    156 |    102 |   34.6% |
| moleditpy\src\moleditpy\modules\mirror_dialog.py        |     70 |      7 |   90.0% |
| moleditpy\src\moleditpy\modules\mol_geometry.py         |    138 |     19 |   86.2% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    273 |     55 |   79.9% |
| moleditpy\src\moleditpy\modules\molecular_scene_handler.py |    846 |    257 |   69.6% |
| moleditpy\src\moleditpy\modules\molecule_scene.py       |    482 |    173 |   64.1% |
| moleditpy\src\moleditpy\modules\move_group_dialog.py    |    376 |    351 |    6.6% |
| moleditpy\src\moleditpy\modules\periodic_table_dialog.py |     33 |     24 |   27.3% |
| moleditpy\src\moleditpy\modules\planarize_dialog.py     |    128 |     37 |   71.1% |
| moleditpy\src\moleditpy\modules\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\plugin_manager.py       |    276 |     68 |   75.4% |
| moleditpy\src\moleditpy\modules\settings_dialog.py      |    713 |    132 |   81.5% |
| moleditpy\src\moleditpy\modules\system_utils.py         |     35 |     18 |   48.6% |
| moleditpy\src\moleditpy\modules\template_preview_item.py |    101 |     77 |   23.8% |
| moleditpy\src\moleditpy\modules\template_preview_view.py |     42 |     14 |   66.7% |
| moleditpy\src\moleditpy\modules\translation_dialog.py   |    196 |    122 |   37.8% |
| moleditpy\src\moleditpy\modules\user_template_dialog.py |    370 |    131 |   64.6% |
| moleditpy\src\moleditpy\modules\zoomable_view.py        |     72 |     39 |   45.8% |
| **TOTAL** | **13934** | **5503** | **60.51%** |

## Test Suite Status
- **Unit tests**: PASSED
- **Integration tests**: PASSED
- **GUI tests**: PASSED

[View Detailed HTML Report](coverage_html/index.html)