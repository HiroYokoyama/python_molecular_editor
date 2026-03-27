# MoleditPy Coverage Report

- **Overall Project Coverage (Full)**: **58.29%**
- **Core Molecular Logic Coverage**: **81.18%**

> [!NOTE]
> **Core Molecular Logic Coverage** excludes UI boilerplate (dialogs, view managers, and interactor styles) to focus on scientific algorithm reliability.

### Core Logic Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      6 |      2 |   66.7% |
| moleditpy\src\moleditpy\core\app_state.py               |    464 |    102 |   78.0% |
| moleditpy\src\moleditpy\core\calculation_worker.py      |    535 |     86 |   83.9% |
| moleditpy\src\moleditpy\core\compute_engine.py          |    482 |    122 |   74.7% |
| moleditpy\src\moleditpy\core\mol_geometry.py            |    238 |     15 |   93.7% |
| moleditpy\src\moleditpy\core\molecular_data.py          |    300 |     63 |   79.0% |
| moleditpy\src\moleditpy\core\molecular_parsers.py       |    349 |     87 |   75.1% |
| moleditpy\src\moleditpy\core\project_io.py              |    232 |     79 |   65.9% |
| moleditpy\src\moleditpy\core\string_importers.py        |    161 |     26 |   83.9% |
| moleditpy\src\moleditpy\plugins\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\plugins\plugin_manager.py       |    272 |     65 |   76.1% |
| moleditpy\src\moleditpy\ui\__init__.py                  |      7 |      3 |   57.1% |
| moleditpy\src\moleditpy\ui\analysis_window.py           |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\ui\main_window.py               |    108 |     19 |   82.4% |
| moleditpy\src\moleditpy\ui\main_window_init.py          |   1158 |    151 |   87.0% |
| moleditpy\src\moleditpy\utils\constants.py              |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\sip_isdeleted_safe.py     |     19 |      8 |   57.9% |
| **TOTAL** | **4533** | **853** | **81.18%** |

### Full Application Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      6 |      2 |   66.7% |
| moleditpy\src\moleditpy\core\app_state.py               |    464 |    102 |   78.0% |
| moleditpy\src\moleditpy\core\calculation_worker.py      |    535 |     86 |   83.9% |
| moleditpy\src\moleditpy\core\compute_engine.py          |    482 |    122 |   74.7% |
| moleditpy\src\moleditpy\core\mol_geometry.py            |    238 |     15 |   93.7% |
| moleditpy\src\moleditpy\core\molecular_data.py          |    300 |     63 |   79.0% |
| moleditpy\src\moleditpy\core\molecular_parsers.py       |    349 |     87 |   75.1% |
| moleditpy\src\moleditpy\core\project_io.py              |    232 |     79 |   65.9% |
| moleditpy\src\moleditpy\core\string_importers.py        |    161 |     26 |   83.9% |
| moleditpy\src\moleditpy\plugins\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\plugins\plugin_manager.py       |    272 |     65 |   76.1% |
| moleditpy\src\moleditpy\ui\__init__.py                  |      7 |      3 |   57.1% |
| moleditpy\src\moleditpy\ui\about_dialog.py              |     64 |     52 |   18.8% |
| moleditpy\src\moleditpy\ui\align_plane_dialog.py        |    162 |     61 |   62.3% |
| moleditpy\src\moleditpy\ui\alignment_dialog.py          |    142 |     47 |   66.9% |
| moleditpy\src\moleditpy\ui\analysis_window.py           |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\ui\angle_dialog.py              |    297 |    144 |   51.5% |
| moleditpy\src\moleditpy\ui\atom_item.py                 |    258 |     34 |   86.8% |
| moleditpy\src\moleditpy\ui\bond_item.py                 |    315 |     78 |   75.2% |
| moleditpy\src\moleditpy\ui\bond_length_dialog.py        |    269 |    136 |   49.4% |
| moleditpy\src\moleditpy\ui\color_settings_dialog.py     |    181 |    170 |    6.1% |
| moleditpy\src\moleditpy\ui\compute_logic.py             |    342 |    148 |   56.7% |
| moleditpy\src\moleditpy\ui\constrained_optimization_dialog.py |    406 |    307 |   24.4% |
| moleditpy\src\moleditpy\ui\custom_interactor_style.py   |    456 |    428 |    6.1% |
| moleditpy\src\moleditpy\ui\custom_qt_interactor.py      |     43 |     35 |   18.6% |
| moleditpy\src\moleditpy\ui\dialog_3d_picking_mixin.py   |     92 |     50 |   45.7% |
| moleditpy\src\moleditpy\ui\dialog_logic.py              |    240 |    157 |   34.6% |
| moleditpy\src\moleditpy\ui\dihedral_dialog.py           |    300 |    139 |   53.7% |
| moleditpy\src\moleditpy\ui\edit_3d.py                   |    218 |    129 |   40.8% |
| moleditpy\src\moleditpy\ui\edit_3d_logic.py             |    235 |    143 |   39.1% |
| moleditpy\src\moleditpy\ui\edit_actions.py              |    845 |    504 |   40.4% |
| moleditpy\src\moleditpy\ui\edit_actions_logic.py        |    773 |    398 |   48.5% |
| moleditpy\src\moleditpy\ui\export_logic.py              |    539 |    162 |   69.9% |
| moleditpy\src\moleditpy\ui\main_window.py               |    108 |     19 |   82.4% |
| moleditpy\src\moleditpy\ui\main_window_init.py          |   1158 |    151 |   87.0% |
| moleditpy\src\moleditpy\ui\mirror_dialog.py             |     70 |      7 |   90.0% |
| moleditpy\src\moleditpy\ui\molecular_scene_handler.py   |    855 |    170 |   80.1% |
| moleditpy\src\moleditpy\ui\molecule_scene.py            |    483 |    173 |   64.2% |
| moleditpy\src\moleditpy\ui\move_group_dialog.py         |    376 |    215 |   42.8% |
| moleditpy\src\moleditpy\ui\periodic_table_dialog.py     |     33 |      7 |   78.8% |
| moleditpy\src\moleditpy\ui\planarize_dialog.py          |    129 |     37 |   71.3% |
| moleditpy\src\moleditpy\ui\settings_dialog.py           |     98 |     38 |   61.2% |
| moleditpy\src\moleditpy\ui\sip_isdeleted_safe.py        |     13 |      4 |   69.2% |
| moleditpy\src\moleditpy\ui\template_preview_item.py     |    101 |     77 |   23.8% |
| moleditpy\src\moleditpy\ui\template_preview_view.py     |     42 |     14 |   66.7% |
| moleditpy\src\moleditpy\ui\translation_dialog.py        |    197 |     87 |   55.8% |
| moleditpy\src\moleditpy\ui\ui_manager.py                |    305 |     84 |   72.5% |
| moleditpy\src\moleditpy\ui\user_template_dialog.py      |    370 |    131 |   64.6% |
| moleditpy\src\moleditpy\ui\view_3d.py                   |    899 |    851 |    5.3% |
| moleditpy\src\moleditpy\ui\view_3d_logic.py             |    916 |    480 |   47.6% |
| moleditpy\src\moleditpy\ui\view_loaders.py              |    155 |    104 |   32.9% |
| moleditpy\src\moleditpy\ui\zoomable_view.py             |     72 |     39 |   45.8% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_2d_tab.py |    100 |      9 |   91.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_3d_tabs.py |    174 |     26 |   85.1% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_other_tab.py |     42 |      1 |   97.6% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_tab_base.py |     11 |      3 |   72.7% |
| moleditpy\src\moleditpy\utils\constants.py              |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\sip_isdeleted_safe.py     |     19 |      8 |   57.9% |
| moleditpy\src\moleditpy\utils\system_utils.py           |     34 |     18 |   47.1% |
| **TOTAL** | **16185** | **6750** | **58.29%** |

## Test Suite Status
- **Unit tests**: PASSED
- **Integration tests**: PASSED
- **GUI tests**: PASSED

[View Detailed HTML Report](coverage_html/index.html)