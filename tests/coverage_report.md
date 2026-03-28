# MoleditPy Coverage Report

- **Overall Project Coverage (Full)**: **40.42%**
- **Core Molecular Logic Coverage**: **75.85%**

> [!NOTE]
> **Core Molecular Logic Coverage** excludes UI boilerplate (dialogs, view managers, and interactor styles) to focus on scientific algorithm reliability.

### Core Logic Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      6 |      2 |   66.7% |
| moleditpy\src\moleditpy\core\mol_geometry.py            |    238 |     17 |   92.9% |
| moleditpy\src\moleditpy\core\molecular_data.py          |    295 |     62 |   79.0% |
| moleditpy\src\moleditpy\plugins\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\plugins\plugin_manager.py       |    272 |     96 |   64.7% |
| moleditpy\src\moleditpy\ui\calculation_worker.py        |    535 |    150 |   72.0% |
| moleditpy\src\moleditpy\ui\compute_engine.py            |    482 |    119 |   75.3% |
| moleditpy\src\moleditpy\ui\molecular_parsers.py         |    349 |     87 |   75.1% |
| moleditpy\src\moleditpy\ui\project_io.py                |    232 |     79 |   65.9% |
| moleditpy\src\moleditpy\ui\string_importers.py          |    161 |     30 |   81.4% |
| moleditpy\src\moleditpy\utils\constants.py              |     31 |      0 |  100.0% |
| **TOTAL** | **2658** | **642** | **75.85%** |

### Full Application Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      6 |      2 |   66.7% |
| moleditpy\src\moleditpy\core\mol_geometry.py            |    238 |     17 |   92.9% |
| moleditpy\src\moleditpy\core\molecular_data.py          |    295 |     62 |   79.0% |
| moleditpy\src\moleditpy\plugins\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\plugins\plugin_manager.py       |    272 |     96 |   64.7% |
| moleditpy\src\moleditpy\ui\__init__.py                  |      7 |      3 |   57.1% |
| moleditpy\src\moleditpy\ui\about_dialog.py              |     64 |     15 |   76.6% |
| moleditpy\src\moleditpy\ui\align_plane_dialog.py        |    162 |     67 |   58.6% |
| moleditpy\src\moleditpy\ui\alignment_dialog.py          |    142 |     51 |   64.1% |
| moleditpy\src\moleditpy\ui\analysis_window.py           |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\ui\angle_dialog.py              |    297 |    168 |   43.4% |
| moleditpy\src\moleditpy\ui\app_state.py                 |    464 |    179 |   61.4% |
| moleditpy\src\moleditpy\ui\atom_item.py                 |    258 |     52 |   79.8% |
| moleditpy\src\moleditpy\ui\bond_item.py                 |    315 |     84 |   73.3% |
| moleditpy\src\moleditpy\ui\bond_length_dialog.py        |    269 |    139 |   48.3% |
| moleditpy\src\moleditpy\ui\calculation_worker.py        |    535 |    150 |   72.0% |
| moleditpy\src\moleditpy\ui\color_settings_dialog.py     |    181 |     53 |   70.7% |
| moleditpy\src\moleditpy\ui\compute_engine.py            |    482 |    119 |   75.3% |
| moleditpy\src\moleditpy\ui\compute_logic.py             |    342 |    302 |   11.7% |
| moleditpy\src\moleditpy\ui\constrained_optimization_dialog.py |    406 |    384 |    5.4% |
| moleditpy\src\moleditpy\ui\custom_interactor_style.py   |    465 |    401 |   13.8% |
| moleditpy\src\moleditpy\ui\custom_qt_interactor.py      |     43 |     35 |   18.6% |
| moleditpy\src\moleditpy\ui\dialog_3d_picking_mixin.py   |     92 |     50 |   45.7% |
| moleditpy\src\moleditpy\ui\dialog_logic.py              |    240 |    197 |   17.9% |
| moleditpy\src\moleditpy\ui\dihedral_dialog.py           |    300 |    139 |   53.7% |
| moleditpy\src\moleditpy\ui\edit_3d.py                   |    218 |    129 |   40.8% |
| moleditpy\src\moleditpy\ui\edit_3d_logic.py             |    235 |    192 |   18.3% |
| moleditpy\src\moleditpy\ui\edit_actions.py              |    845 |    504 |   40.4% |
| moleditpy\src\moleditpy\ui\edit_actions_logic.py        |    773 |    716 |    7.4% |
| moleditpy\src\moleditpy\ui\export_logic.py              |    539 |    162 |   69.9% |
| moleditpy\src\moleditpy\ui\main_window.py               |    108 |     63 |   41.7% |
| moleditpy\src\moleditpy\ui\main_window_init.py          |   1158 |   1078 |    6.9% |
| moleditpy\src\moleditpy\ui\mirror_dialog.py             |     70 |      8 |   88.6% |
| moleditpy\src\moleditpy\ui\molecular_parsers.py         |    349 |     87 |   75.1% |
| moleditpy\src\moleditpy\ui\molecular_scene_handler.py   |    855 |    281 |   67.1% |
| moleditpy\src\moleditpy\ui\molecule_scene.py            |    483 |    189 |   60.9% |
| moleditpy\src\moleditpy\ui\move_group_dialog.py         |    376 |    299 |   20.5% |
| moleditpy\src\moleditpy\ui\periodic_table_dialog.py     |     33 |     24 |   27.3% |
| moleditpy\src\moleditpy\ui\planarize_dialog.py          |    129 |     37 |   71.3% |
| moleditpy\src\moleditpy\ui\project_io.py                |    232 |     79 |   65.9% |
| moleditpy\src\moleditpy\ui\settings_dialog.py           |     98 |     82 |   16.3% |
| moleditpy\src\moleditpy\ui\sip_isdeleted_safe.py        |     13 |      4 |   69.2% |
| moleditpy\src\moleditpy\ui\string_importers.py          |    161 |     30 |   81.4% |
| moleditpy\src\moleditpy\ui\template_preview_item.py     |    101 |     77 |   23.8% |
| moleditpy\src\moleditpy\ui\template_preview_view.py     |     42 |     33 |   21.4% |
| moleditpy\src\moleditpy\ui\translation_dialog.py        |    197 |    132 |   33.0% |
| moleditpy\src\moleditpy\ui\ui_manager.py                |    305 |    225 |   26.2% |
| moleditpy\src\moleditpy\ui\user_template_dialog.py      |    370 |    339 |    8.4% |
| moleditpy\src\moleditpy\ui\view_3d.py                   |    899 |    703 |   21.8% |
| moleditpy\src\moleditpy\ui\view_3d_logic.py             |    916 |    866 |    5.5% |
| moleditpy\src\moleditpy\ui\view_loaders.py              |    155 |    140 |    9.7% |
| moleditpy\src\moleditpy\ui\zoomable_view.py             |     72 |     63 |   12.5% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_2d_tab.py |    100 |     85 |   15.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_3d_tabs.py |    174 |    154 |   11.5% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_other_tab.py |     42 |     34 |   19.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_tab_base.py |     11 |      5 |   54.5% |
| moleditpy\src\moleditpy\utils\constants.py              |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\sip_isdeleted_safe.py     |     19 |      8 |   57.9% |
| moleditpy\src\moleditpy\utils\system_utils.py           |     34 |     27 |   20.6% |
| **TOTAL** | **16189** | **9645** | **40.42%** |

## Test Suite Status
- **Unit tests**: PASSED
- **Integration tests**: PASSED
- **GUI tests**: PASSED

[View Detailed HTML Report](coverage_html/index.html)