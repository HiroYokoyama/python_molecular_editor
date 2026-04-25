# MoleditPy Coverage Report

- **Overall Project Coverage (Full)**: **74.47%**
- **Core Molecular Logic Coverage**: **80.50%**

> [!NOTE]
> **Core Molecular Logic Coverage** excludes UI boilerplate (dialogs, view managers, and interactor styles) to focus on scientific algorithm reliability.

### Core Logic Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      6 |      2 |   66.7% |
| moleditpy\src\moleditpy\main.py                         |     79 |     22 |   72.2% |
| moleditpy\src\moleditpy\core\mol_geometry.py            |    260 |     19 |   92.7% |
| moleditpy\src\moleditpy\core\molecular_data.py          |    307 |     38 |   87.6% |
| moleditpy\src\moleditpy\plugins\plugin_interface.py     |    127 |      0 |  100.0% |
| moleditpy\src\moleditpy\plugins\plugin_manager.py       |    370 |     15 |   95.9% |
| moleditpy\src\moleditpy\plugins\plugin_manager_window.py |    179 |      0 |  100.0% |
| moleditpy\src\moleditpy\ui\about_dialog.py              |     65 |     15 |   76.9% |
| moleditpy\src\moleditpy\ui\align_plane_dialog.py        |    123 |     17 |   86.2% |
| moleditpy\src\moleditpy\ui\alignment_dialog.py          |    140 |     13 |   90.7% |
| moleditpy\src\moleditpy\ui\angle_dialog.py              |    228 |     18 |   92.1% |
| moleditpy\src\moleditpy\ui\app_state.py                 |    422 |    102 |   75.8% |
| moleditpy\src\moleditpy\ui\atom_item.py                 |    265 |     34 |   87.2% |
| moleditpy\src\moleditpy\ui\bond_item.py                 |    301 |     71 |   76.4% |
| moleditpy\src\moleditpy\ui\bond_length_dialog.py        |    207 |     13 |   93.7% |
| moleditpy\src\moleditpy\ui\calculation_worker.py        |    547 |     91 |   83.4% |
| moleditpy\src\moleditpy\ui\color_settings_dialog.py     |    197 |     65 |   67.0% |
| moleditpy\src\moleditpy\ui\compute_logic.py             |    401 |     96 |   76.1% |
| moleditpy\src\moleditpy\ui\dialog_logic.py              |    228 |     31 |   86.4% |
| moleditpy\src\moleditpy\ui\dihedral_dialog.py           |    223 |     17 |   92.4% |
| moleditpy\src\moleditpy\ui\edit_3d_logic.py             |    238 |     88 |   63.0% |
| moleditpy\src\moleditpy\ui\edit_actions_logic.py        |    852 |    215 |   74.8% |
| moleditpy\src\moleditpy\ui\export_logic.py              |    512 |    150 |   70.7% |
| moleditpy\src\moleditpy\ui\io_logic.py                  |    637 |    131 |   79.4% |
| moleditpy\src\moleditpy\ui\main_window_init.py          |   1190 |    167 |   86.0% |
| moleditpy\src\moleditpy\ui\molecular_scene_handler.py   |    858 |    178 |   79.3% |
| moleditpy\src\moleditpy\ui\molecule_scene.py            |    502 |    147 |   70.7% |
| moleditpy\src\moleditpy\ui\string_importers.py          |    128 |     27 |   78.9% |
| moleditpy\src\moleditpy\ui\ui_manager.py                |    373 |    112 |   70.0% |
| moleditpy\src\moleditpy\ui\view_3d_logic.py             |    992 |    246 |   75.2% |
| moleditpy\src\moleditpy\utils\constants.py              |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\default_settings.py       |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\sip_isdeleted_safe.py     |     19 |      7 |   63.2% |
| **TOTAL** | **11008** | **2147** | **80.50%** |

### Full Application Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      6 |      2 |   66.7% |
| moleditpy\src\moleditpy\main.py                         |     79 |     22 |   72.2% |
| moleditpy\src\moleditpy\core\mol_geometry.py            |    260 |     19 |   92.7% |
| moleditpy\src\moleditpy\core\molecular_data.py          |    307 |     38 |   87.6% |
| moleditpy\src\moleditpy\plugins\plugin_interface.py     |    127 |      0 |  100.0% |
| moleditpy\src\moleditpy\plugins\plugin_manager.py       |    370 |     15 |   95.9% |
| moleditpy\src\moleditpy\plugins\plugin_manager_window.py |    179 |      0 |  100.0% |
| moleditpy\src\moleditpy\ui\__init__.py                  |      7 |      3 |   57.1% |
| moleditpy\src\moleditpy\ui\about_dialog.py              |     65 |     15 |   76.9% |
| moleditpy\src\moleditpy\ui\align_plane_dialog.py        |    125 |     18 |   85.6% |
| moleditpy\src\moleditpy\ui\alignment_dialog.py          |    142 |     14 |   90.1% |
| moleditpy\src\moleditpy\ui\analysis_window.py           |    116 |     25 |   78.4% |
| moleditpy\src\moleditpy\ui\angle_dialog.py              |    230 |     19 |   91.7% |
| moleditpy\src\moleditpy\ui\app_state.py                 |    422 |    102 |   75.8% |
| moleditpy\src\moleditpy\ui\atom_item.py                 |    265 |     34 |   87.2% |
| moleditpy\src\moleditpy\ui\base_picking_dialog.py       |     68 |     20 |   70.6% |
| moleditpy\src\moleditpy\ui\bond_item.py                 |    301 |     71 |   76.4% |
| moleditpy\src\moleditpy\ui\bond_length_dialog.py        |    209 |     14 |   93.3% |
| moleditpy\src\moleditpy\ui\calculation_worker.py        |    547 |     91 |   83.4% |
| moleditpy\src\moleditpy\ui\color_settings_dialog.py     |    197 |     65 |   67.0% |
| moleditpy\src\moleditpy\ui\compute_logic.py             |    401 |     96 |   76.1% |
| moleditpy\src\moleditpy\ui\constrained_optimization_dialog.py |    408 |    121 |   70.3% |
| moleditpy\src\moleditpy\ui\custom_interactor_style.py   |    454 |    391 |   13.9% |
| moleditpy\src\moleditpy\ui\custom_qt_interactor.py      |     43 |     35 |   18.6% |
| moleditpy\src\moleditpy\ui\dialog_3d_picking_mixin.py   |    117 |     60 |   48.7% |
| moleditpy\src\moleditpy\ui\dialog_logic.py              |    228 |     31 |   86.4% |
| moleditpy\src\moleditpy\ui\dihedral_dialog.py           |    225 |     18 |   92.0% |
| moleditpy\src\moleditpy\ui\edit_3d_logic.py             |    238 |     88 |   63.0% |
| moleditpy\src\moleditpy\ui\edit_actions_logic.py        |    852 |    215 |   74.8% |
| moleditpy\src\moleditpy\ui\export_logic.py              |    512 |    150 |   70.7% |
| moleditpy\src\moleditpy\ui\geometry_base_dialog.py      |     61 |      8 |   86.9% |
| moleditpy\src\moleditpy\ui\io_logic.py                  |    637 |    131 |   79.4% |
| moleditpy\src\moleditpy\ui\main_window.py               |     66 |     20 |   69.7% |
| moleditpy\src\moleditpy\ui\main_window_init.py          |   1190 |    167 |   86.0% |
| moleditpy\src\moleditpy\ui\mirror_dialog.py             |     70 |      7 |   90.0% |
| moleditpy\src\moleditpy\ui\molecular_scene_handler.py   |    858 |    178 |   79.3% |
| moleditpy\src\moleditpy\ui\molecule_scene.py            |    502 |    147 |   70.7% |
| moleditpy\src\moleditpy\ui\move_group_dialog.py         |    352 |    156 |   55.7% |
| moleditpy\src\moleditpy\ui\periodic_table_dialog.py     |     33 |      7 |   78.8% |
| moleditpy\src\moleditpy\ui\planarize_dialog.py          |    106 |     11 |   89.6% |
| moleditpy\src\moleditpy\ui\settings_dialog.py           |    107 |     89 |   16.8% |
| moleditpy\src\moleditpy\ui\string_importers.py          |    128 |     27 |   78.9% |
| moleditpy\src\moleditpy\ui\template_preview_item.py     |    103 |     77 |   25.2% |
| moleditpy\src\moleditpy\ui\template_preview_view.py     |     45 |     14 |   68.9% |
| moleditpy\src\moleditpy\ui\translation_dialog.py        |    237 |      9 |   96.2% |
| moleditpy\src\moleditpy\ui\ui_manager.py                |    373 |    112 |   70.0% |
| moleditpy\src\moleditpy\ui\user_template_dialog.py      |    365 |    128 |   64.9% |
| moleditpy\src\moleditpy\ui\view_3d_logic.py             |    992 |    246 |   75.2% |
| moleditpy\src\moleditpy\ui\zoomable_view.py             |     73 |     37 |   49.3% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_2d_tab.py |     81 |     68 |   16.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_3d_tabs.py |    135 |    117 |   13.3% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_other_tab.py |     54 |     42 |   22.2% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_tab_base.py |     36 |     24 |   33.3% |
| moleditpy\src\moleditpy\utils\constants.py              |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\default_settings.py       |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\sip_isdeleted_safe.py     |     19 |      7 |   63.2% |
| moleditpy\src\moleditpy\utils\system_utils.py           |     36 |      2 |   94.4% |
| **TOTAL** | **14191** | **3623** | **74.47%** |

## Test Suite Status
- **Unit tests**: PASSED
- **Integration tests**: PASSED
- **GUI tests**: PASSED

[View Detailed HTML Report](coverage_html/index.html)