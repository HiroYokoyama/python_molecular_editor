# MoleditPy Coverage Report

- **Overall Project Coverage (Full)**: **58.39%**
- **Core Molecular Logic Coverage**: **81.19%**

> [!NOTE]
> **Core Molecular Logic Coverage** excludes UI boilerplate (dialogs, view managers, and interactor styles) to focus on scientific algorithm reliability.

### Core Logic Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      6 |      2 |   66.7% |
| moleditpy\src\moleditpy\core\app_state.py               |    465 |    102 |   78.1% |
| moleditpy\src\moleditpy\core\calculation_worker.py      |    535 |     86 |   83.9% |
| moleditpy\src\moleditpy\core\compute_engine.py          |    483 |    122 |   74.7% |
| moleditpy\src\moleditpy\core\mol_geometry.py            |    239 |     15 |   93.7% |
| moleditpy\src\moleditpy\core\molecular_data.py          |    301 |     63 |   79.1% |
| moleditpy\src\moleditpy\core\molecular_parsers.py       |    350 |     87 |   75.1% |
| moleditpy\src\moleditpy\core\project_io.py              |    234 |     79 |   66.2% |
| moleditpy\src\moleditpy\core\string_importers.py        |    162 |     26 |   84.0% |
| moleditpy\src\moleditpy\plugins\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\plugins\plugin_manager.py       |    273 |     65 |   76.2% |
| moleditpy\src\moleditpy\ui\__init__.py                  |      7 |      3 |   57.1% |
| moleditpy\src\moleditpy\ui\analysis_window.py           |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\ui\main_window.py               |    110 |     21 |   80.9% |
| moleditpy\src\moleditpy\ui\main_window_init.py          |   1160 |    151 |   87.0% |
| moleditpy\src\moleditpy\utils\constants.py              |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\sip_isdeleted_safe.py     |     19 |      8 |   57.9% |
| **TOTAL** | **4546** | **855** | **81.19%** |

### Full Application Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      6 |      2 |   66.7% |
| moleditpy\src\moleditpy\core\app_state.py               |    465 |    102 |   78.1% |
| moleditpy\src\moleditpy\core\calculation_worker.py      |    535 |     86 |   83.9% |
| moleditpy\src\moleditpy\core\compute_engine.py          |    483 |    122 |   74.7% |
| moleditpy\src\moleditpy\core\mol_geometry.py            |    239 |     15 |   93.7% |
| moleditpy\src\moleditpy\core\molecular_data.py          |    301 |     63 |   79.1% |
| moleditpy\src\moleditpy\core\molecular_parsers.py       |    350 |     87 |   75.1% |
| moleditpy\src\moleditpy\core\project_io.py              |    234 |     79 |   66.2% |
| moleditpy\src\moleditpy\core\string_importers.py        |    162 |     26 |   84.0% |
| moleditpy\src\moleditpy\plugins\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\plugins\plugin_manager.py       |    273 |     65 |   76.2% |
| moleditpy\src\moleditpy\ui\__init__.py                  |      7 |      3 |   57.1% |
| moleditpy\src\moleditpy\ui\about_dialog.py              |     65 |     52 |   20.0% |
| moleditpy\src\moleditpy\ui\align_plane_dialog.py        |    162 |     61 |   62.3% |
| moleditpy\src\moleditpy\ui\alignment_dialog.py          |    142 |     47 |   66.9% |
| moleditpy\src\moleditpy\ui\analysis_window.py           |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\ui\angle_dialog.py              |    298 |    144 |   51.7% |
| moleditpy\src\moleditpy\ui\atom_item.py                 |    258 |     34 |   86.8% |
| moleditpy\src\moleditpy\ui\bond_item.py                 |    311 |     74 |   76.2% |
| moleditpy\src\moleditpy\ui\bond_length_dialog.py        |    270 |    136 |   49.6% |
| moleditpy\src\moleditpy\ui\color_settings_dialog.py     |    181 |    170 |    6.1% |
| moleditpy\src\moleditpy\ui\compute_logic.py             |    344 |    148 |   57.0% |
| moleditpy\src\moleditpy\ui\constrained_optimization_dialog.py |    407 |    307 |   24.6% |
| moleditpy\src\moleditpy\ui\custom_interactor_style.py   |    457 |    428 |    6.3% |
| moleditpy\src\moleditpy\ui\custom_qt_interactor.py      |     43 |     35 |   18.6% |
| moleditpy\src\moleditpy\ui\dialog_3d_picking_mixin.py   |     92 |     50 |   45.7% |
| moleditpy\src\moleditpy\ui\dialog_logic.py              |    243 |    157 |   35.4% |
| moleditpy\src\moleditpy\ui\dihedral_dialog.py           |    301 |    139 |   53.8% |
| moleditpy\src\moleditpy\ui\edit_3d.py                   |    219 |    129 |   41.1% |
| moleditpy\src\moleditpy\ui\edit_3d_logic.py             |    235 |    143 |   39.1% |
| moleditpy\src\moleditpy\ui\edit_actions.py              |    846 |    504 |   40.4% |
| moleditpy\src\moleditpy\ui\edit_actions_logic.py        |    777 |    398 |   48.8% |
| moleditpy\src\moleditpy\ui\export_logic.py              |    540 |    162 |   70.0% |
| moleditpy\src\moleditpy\ui\main_window.py               |    110 |     21 |   80.9% |
| moleditpy\src\moleditpy\ui\main_window_init.py          |   1160 |    151 |   87.0% |
| moleditpy\src\moleditpy\ui\mirror_dialog.py             |     70 |      7 |   90.0% |
| moleditpy\src\moleditpy\ui\molecular_scene_handler.py   |    855 |    170 |   80.1% |
| moleditpy\src\moleditpy\ui\molecule_scene.py            |    483 |    173 |   64.2% |
| moleditpy\src\moleditpy\ui\move_group_dialog.py         |    377 |    215 |   43.0% |
| moleditpy\src\moleditpy\ui\periodic_table_dialog.py     |     33 |      7 |   78.8% |
| moleditpy\src\moleditpy\ui\planarize_dialog.py          |    129 |     37 |   71.3% |
| moleditpy\src\moleditpy\ui\settings_dialog.py           |     97 |     36 |   62.9% |
| moleditpy\src\moleditpy\ui\sip_isdeleted_safe.py        |     13 |      4 |   69.2% |
| moleditpy\src\moleditpy\ui\template_preview_item.py     |    101 |     77 |   23.8% |
| moleditpy\src\moleditpy\ui\template_preview_view.py     |     42 |     14 |   66.7% |
| moleditpy\src\moleditpy\ui\translation_dialog.py        |    197 |     87 |   55.8% |
| moleditpy\src\moleditpy\ui\ui_manager.py                |    306 |     84 |   72.5% |
| moleditpy\src\moleditpy\ui\user_template_dialog.py      |    370 |    131 |   64.6% |
| moleditpy\src\moleditpy\ui\view_3d.py                   |    900 |    851 |    5.4% |
| moleditpy\src\moleditpy\ui\view_3d_logic.py             |    917 |    480 |   47.7% |
| moleditpy\src\moleditpy\ui\view_loaders.py              |    156 |    104 |   33.3% |
| moleditpy\src\moleditpy\ui\zoomable_view.py             |     72 |     39 |   45.8% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_2d_tab.py |     99 |      9 |   90.9% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_3d_tabs.py |    173 |     26 |   85.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_other_tab.py |     42 |      1 |   97.6% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_tab_base.py |     11 |      3 |   72.7% |
| moleditpy\src\moleditpy\utils\constants.py              |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\sip_isdeleted_safe.py     |     19 |      8 |   57.9% |
| moleditpy\src\moleditpy\utils\system_utils.py           |     34 |     18 |   47.1% |
| **TOTAL** | **16214** | **6746** | **58.39%** |

## Test Suite Status
- **Unit tests**: PASSED
- **Integration tests**: PASSED
- **GUI tests**: PASSED

[View Detailed HTML Report](coverage_html/index.html)