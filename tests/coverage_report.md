# MoleditPy Coverage Report

- **Overall Project Coverage (Full)**: **74.42%**
- **Core Molecular Logic Coverage**: **80.49%**

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
| moleditpy\src\moleditpy\plugins\plugin_manager_window.py |    178 |      0 |  100.0% |
| moleditpy\src\moleditpy\ui\about_dialog.py              |     65 |     15 |   76.9% |
| moleditpy\src\moleditpy\ui\align_plane_dialog.py        |    121 |     17 |   86.0% |
| moleditpy\src\moleditpy\ui\alignment_dialog.py          |    138 |     13 |   90.6% |
| moleditpy\src\moleditpy\ui\angle_dialog.py              |    226 |     18 |   92.0% |
| moleditpy\src\moleditpy\ui\app_state.py                 |    422 |    102 |   75.8% |
| moleditpy\src\moleditpy\ui\atom_item.py                 |    265 |     34 |   87.2% |
| moleditpy\src\moleditpy\ui\bond_item.py                 |    301 |     71 |   76.4% |
| moleditpy\src\moleditpy\ui\bond_length_dialog.py        |    205 |     13 |   93.7% |
| moleditpy\src\moleditpy\ui\calculation_worker.py        |    546 |     91 |   83.3% |
| moleditpy\src\moleditpy\ui\color_settings_dialog.py     |    197 |     65 |   67.0% |
| moleditpy\src\moleditpy\ui\compute_logic.py             |    401 |     96 |   76.1% |
| moleditpy\src\moleditpy\ui\dialog_logic.py              |    228 |     31 |   86.4% |
| moleditpy\src\moleditpy\ui\dihedral_dialog.py           |    221 |     17 |   92.3% |
| moleditpy\src\moleditpy\ui\edit_3d_logic.py             |    238 |     88 |   63.0% |
| moleditpy\src\moleditpy\ui\edit_actions_logic.py        |    852 |    215 |   74.8% |
| moleditpy\src\moleditpy\ui\export_logic.py              |    512 |    150 |   70.7% |
| moleditpy\src\moleditpy\ui\io_logic.py                  |    637 |    131 |   79.4% |
| moleditpy\src\moleditpy\ui\main_window_init.py          |   1190 |    167 |   86.0% |
| moleditpy\src\moleditpy\ui\molecular_scene_handler.py   |    858 |    178 |   79.3% |
| moleditpy\src\moleditpy\ui\molecule_scene.py            |    502 |    145 |   71.1% |
| moleditpy\src\moleditpy\ui\string_importers.py          |    128 |     27 |   78.9% |
| moleditpy\src\moleditpy\ui\ui_manager.py                |    372 |    112 |   69.9% |
| moleditpy\src\moleditpy\ui\view_3d_logic.py             |    992 |    246 |   75.2% |
| moleditpy\src\moleditpy\utils\constants.py              |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\default_settings.py       |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\sip_isdeleted_safe.py     |     19 |      7 |   63.2% |
| **TOTAL** | **10995** | **2145** | **80.49%** |

### Full Application Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      6 |      2 |   66.7% |
| moleditpy\src\moleditpy\main.py                         |     79 |     22 |   72.2% |
| moleditpy\src\moleditpy\core\mol_geometry.py            |    260 |     19 |   92.7% |
| moleditpy\src\moleditpy\core\molecular_data.py          |    307 |     38 |   87.6% |
| moleditpy\src\moleditpy\plugins\plugin_interface.py     |    127 |      0 |  100.0% |
| moleditpy\src\moleditpy\plugins\plugin_manager.py       |    370 |     15 |   95.9% |
| moleditpy\src\moleditpy\plugins\plugin_manager_window.py |    178 |      0 |  100.0% |
| moleditpy\src\moleditpy\ui\__init__.py                  |      7 |      3 |   57.1% |
| moleditpy\src\moleditpy\ui\about_dialog.py              |     65 |     15 |   76.9% |
| moleditpy\src\moleditpy\ui\align_plane_dialog.py        |    121 |     17 |   86.0% |
| moleditpy\src\moleditpy\ui\alignment_dialog.py          |    138 |     13 |   90.6% |
| moleditpy\src\moleditpy\ui\analysis_window.py           |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\ui\angle_dialog.py              |    226 |     18 |   92.0% |
| moleditpy\src\moleditpy\ui\app_state.py                 |    422 |    102 |   75.8% |
| moleditpy\src\moleditpy\ui\atom_item.py                 |    265 |     34 |   87.2% |
| moleditpy\src\moleditpy\ui\base_picking_dialog.py       |     63 |     19 |   69.8% |
| moleditpy\src\moleditpy\ui\bond_item.py                 |    301 |     71 |   76.4% |
| moleditpy\src\moleditpy\ui\bond_length_dialog.py        |    205 |     13 |   93.7% |
| moleditpy\src\moleditpy\ui\calculation_worker.py        |    546 |     91 |   83.3% |
| moleditpy\src\moleditpy\ui\color_settings_dialog.py     |    197 |     65 |   67.0% |
| moleditpy\src\moleditpy\ui\compute_logic.py             |    401 |     96 |   76.1% |
| moleditpy\src\moleditpy\ui\constrained_optimization_dialog.py |    406 |    121 |   70.2% |
| moleditpy\src\moleditpy\ui\custom_interactor_style.py   |    454 |    391 |   13.9% |
| moleditpy\src\moleditpy\ui\custom_qt_interactor.py      |     43 |     35 |   18.6% |
| moleditpy\src\moleditpy\ui\dialog_3d_picking_mixin.py   |    115 |     60 |   47.8% |
| moleditpy\src\moleditpy\ui\dialog_logic.py              |    228 |     31 |   86.4% |
| moleditpy\src\moleditpy\ui\dihedral_dialog.py           |    221 |     17 |   92.3% |
| moleditpy\src\moleditpy\ui\edit_3d_logic.py             |    238 |     88 |   63.0% |
| moleditpy\src\moleditpy\ui\edit_actions_logic.py        |    852 |    215 |   74.8% |
| moleditpy\src\moleditpy\ui\export_logic.py              |    512 |    150 |   70.7% |
| moleditpy\src\moleditpy\ui\geometry_base_dialog.py      |     56 |      7 |   87.5% |
| moleditpy\src\moleditpy\ui\io_logic.py                  |    637 |    131 |   79.4% |
| moleditpy\src\moleditpy\ui\main_window.py               |     66 |     20 |   69.7% |
| moleditpy\src\moleditpy\ui\main_window_init.py          |   1190 |    167 |   86.0% |
| moleditpy\src\moleditpy\ui\mirror_dialog.py             |     70 |      7 |   90.0% |
| moleditpy\src\moleditpy\ui\molecular_scene_handler.py   |    858 |    178 |   79.3% |
| moleditpy\src\moleditpy\ui\molecule_scene.py            |    502 |    145 |   71.1% |
| moleditpy\src\moleditpy\ui\move_group_dialog.py         |    350 |    156 |   55.4% |
| moleditpy\src\moleditpy\ui\periodic_table_dialog.py     |     33 |      7 |   78.8% |
| moleditpy\src\moleditpy\ui\planarize_dialog.py          |    102 |     10 |   90.2% |
| moleditpy\src\moleditpy\ui\settings_dialog.py           |    106 |     89 |   16.0% |
| moleditpy\src\moleditpy\ui\string_importers.py          |    128 |     27 |   78.9% |
| moleditpy\src\moleditpy\ui\template_preview_item.py     |    101 |     77 |   23.8% |
| moleditpy\src\moleditpy\ui\template_preview_view.py     |     42 |     14 |   66.7% |
| moleditpy\src\moleditpy\ui\translation_dialog.py        |    236 |      9 |   96.2% |
| moleditpy\src\moleditpy\ui\ui_manager.py                |    372 |    112 |   69.9% |
| moleditpy\src\moleditpy\ui\user_template_dialog.py      |    364 |    128 |   64.8% |
| moleditpy\src\moleditpy\ui\view_3d_logic.py             |    992 |    246 |   75.2% |
| moleditpy\src\moleditpy\ui\zoomable_view.py             |     70 |     37 |   47.1% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_2d_tab.py |     79 |     68 |   13.9% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_3d_tabs.py |    133 |    117 |   12.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_other_tab.py |     52 |     42 |   19.2% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_tab_base.py |     34 |     24 |   29.4% |
| moleditpy\src\moleditpy\utils\constants.py              |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\default_settings.py       |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\sip_isdeleted_safe.py     |     19 |      7 |   63.2% |
| moleditpy\src\moleditpy\utils\system_utils.py           |     36 |      2 |   94.4% |
| **TOTAL** | **14127** | **3613** | **74.42%** |

## Test Suite Status
- **Unit tests**: PASSED
- **Integration tests**: PASSED
- **GUI tests**: PASSED

[View Detailed HTML Report](coverage_html/index.html)