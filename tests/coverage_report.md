# MoleditPy Coverage Report

- **Overall Project Coverage (Full)**: **67.41%**
- **Core Molecular Logic Coverage**: **75.30%**

> [!NOTE]
> **Core Molecular Logic Coverage** excludes UI boilerplate (dialogs, view managers, and interactor styles) to focus on scientific algorithm reliability.

### Core Logic Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      6 |      2 |   66.7% |
| moleditpy\src\moleditpy\core\mol_geometry.py            |    260 |     19 |   92.7% |
| moleditpy\src\moleditpy\core\molecular_data.py          |    302 |     55 |   81.8% |
| moleditpy\src\moleditpy\plugins\plugin_interface.py     |    118 |     20 |   83.1% |
| moleditpy\src\moleditpy\plugins\plugin_manager.py       |    325 |    108 |   66.8% |
| moleditpy\src\moleditpy\ui\about_dialog.py              |     64 |     15 |   76.6% |
| moleditpy\src\moleditpy\ui\align_plane_dialog.py        |    122 |     42 |   65.6% |
| moleditpy\src\moleditpy\ui\alignment_dialog.py          |    142 |     47 |   66.9% |
| moleditpy\src\moleditpy\ui\angle_dialog.py              |    227 |     94 |   58.6% |
| moleditpy\src\moleditpy\ui\app_state.py                 |    421 |    102 |   75.8% |
| moleditpy\src\moleditpy\ui\atom_item.py                 |    265 |     34 |   87.2% |
| moleditpy\src\moleditpy\ui\bond_item.py                 |    306 |     72 |   76.5% |
| moleditpy\src\moleditpy\ui\bond_length_dialog.py        |    206 |     87 |   57.8% |
| moleditpy\src\moleditpy\ui\calculation_worker.py        |    546 |     91 |   83.3% |
| moleditpy\src\moleditpy\ui\color_settings_dialog.py     |    197 |     65 |   67.0% |
| moleditpy\src\moleditpy\ui\compute_logic.py             |    401 |     96 |   76.1% |
| moleditpy\src\moleditpy\ui\dialog_logic.py              |    232 |     33 |   85.8% |
| moleditpy\src\moleditpy\ui\dihedral_dialog.py           |    222 |     90 |   59.5% |
| moleditpy\src\moleditpy\ui\edit_3d_logic.py             |    236 |    128 |   45.8% |
| moleditpy\src\moleditpy\ui\edit_actions_logic.py        |    852 |    215 |   74.8% |
| moleditpy\src\moleditpy\ui\export_logic.py              |    512 |    150 |   70.7% |
| moleditpy\src\moleditpy\ui\io_logic.py                  |    643 |    131 |   79.6% |
| moleditpy\src\moleditpy\ui\main_window_init.py          |   1186 |    165 |   86.1% |
| moleditpy\src\moleditpy\ui\molecular_scene_handler.py   |    860 |    185 |   78.5% |
| moleditpy\src\moleditpy\ui\molecule_scene.py            |    506 |    187 |   63.0% |
| moleditpy\src\moleditpy\ui\sip_isdeleted_safe.py        |     13 |      4 |   69.2% |
| moleditpy\src\moleditpy\ui\string_importers.py          |    165 |     26 |   84.2% |
| moleditpy\src\moleditpy\ui\ui_manager.py                |    372 |    137 |   63.2% |
| moleditpy\src\moleditpy\ui\view_3d_logic.py             |    986 |    246 |   75.1% |
| moleditpy\src\moleditpy\utils\constants.py              |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\default_settings.py       |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\sip_isdeleted_safe.py     |     19 |      8 |   57.9% |
| **TOTAL** | **10744** | **2654** | **75.30%** |

### Full Application Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      6 |      2 |   66.7% |
| moleditpy\src\moleditpy\core\mol_geometry.py            |    260 |     19 |   92.7% |
| moleditpy\src\moleditpy\core\molecular_data.py          |    302 |     55 |   81.8% |
| moleditpy\src\moleditpy\plugins\plugin_interface.py     |    118 |     20 |   83.1% |
| moleditpy\src\moleditpy\plugins\plugin_manager.py       |    325 |    108 |   66.8% |
| moleditpy\src\moleditpy\ui\__init__.py                  |      7 |      3 |   57.1% |
| moleditpy\src\moleditpy\ui\about_dialog.py              |     64 |     15 |   76.6% |
| moleditpy\src\moleditpy\ui\align_plane_dialog.py        |    122 |     42 |   65.6% |
| moleditpy\src\moleditpy\ui\alignment_dialog.py          |    142 |     47 |   66.9% |
| moleditpy\src\moleditpy\ui\analysis_window.py           |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\ui\angle_dialog.py              |    227 |     94 |   58.6% |
| moleditpy\src\moleditpy\ui\app_state.py                 |    421 |    102 |   75.8% |
| moleditpy\src\moleditpy\ui\atom_item.py                 |    265 |     34 |   87.2% |
| moleditpy\src\moleditpy\ui\base_picking_dialog.py       |     63 |     19 |   69.8% |
| moleditpy\src\moleditpy\ui\bond_item.py                 |    306 |     72 |   76.5% |
| moleditpy\src\moleditpy\ui\bond_length_dialog.py        |    206 |     87 |   57.8% |
| moleditpy\src\moleditpy\ui\calculation_worker.py        |    546 |     91 |   83.3% |
| moleditpy\src\moleditpy\ui\color_settings_dialog.py     |    197 |     65 |   67.0% |
| moleditpy\src\moleditpy\ui\compute_logic.py             |    401 |     96 |   76.1% |
| moleditpy\src\moleditpy\ui\constrained_optimization_dialog.py |    406 |    307 |   24.4% |
| moleditpy\src\moleditpy\ui\custom_interactor_style.py   |    456 |    392 |   14.0% |
| moleditpy\src\moleditpy\ui\custom_qt_interactor.py      |     43 |     35 |   18.6% |
| moleditpy\src\moleditpy\ui\dialog_3d_picking_mixin.py   |     95 |     52 |   45.3% |
| moleditpy\src\moleditpy\ui\dialog_logic.py              |    232 |     33 |   85.8% |
| moleditpy\src\moleditpy\ui\dihedral_dialog.py           |    222 |     90 |   59.5% |
| moleditpy\src\moleditpy\ui\edit_3d_logic.py             |    236 |    128 |   45.8% |
| moleditpy\src\moleditpy\ui\edit_actions_logic.py        |    852 |    215 |   74.8% |
| moleditpy\src\moleditpy\ui\export_logic.py              |    512 |    150 |   70.7% |
| moleditpy\src\moleditpy\ui\geometry_base_dialog.py      |     56 |     41 |   26.8% |
| moleditpy\src\moleditpy\ui\io_logic.py                  |    643 |    131 |   79.6% |
| moleditpy\src\moleditpy\ui\main_window.py               |     64 |     20 |   68.8% |
| moleditpy\src\moleditpy\ui\main_window_init.py          |   1186 |    165 |   86.1% |
| moleditpy\src\moleditpy\ui\mirror_dialog.py             |     70 |      7 |   90.0% |
| moleditpy\src\moleditpy\ui\molecular_scene_handler.py   |    860 |    185 |   78.5% |
| moleditpy\src\moleditpy\ui\molecule_scene.py            |    506 |    187 |   63.0% |
| moleditpy\src\moleditpy\ui\move_group_dialog.py         |    350 |    208 |   40.6% |
| moleditpy\src\moleditpy\ui\periodic_table_dialog.py     |     33 |      7 |   78.8% |
| moleditpy\src\moleditpy\ui\planarize_dialog.py          |    103 |     28 |   72.8% |
| moleditpy\src\moleditpy\ui\settings_dialog.py           |    107 |     89 |   16.8% |
| moleditpy\src\moleditpy\ui\sip_isdeleted_safe.py        |     13 |      4 |   69.2% |
| moleditpy\src\moleditpy\ui\string_importers.py          |    165 |     26 |   84.2% |
| moleditpy\src\moleditpy\ui\template_preview_item.py     |    101 |     77 |   23.8% |
| moleditpy\src\moleditpy\ui\template_preview_view.py     |     42 |     14 |   66.7% |
| moleditpy\src\moleditpy\ui\translation_dialog.py        |    108 |     35 |   67.6% |
| moleditpy\src\moleditpy\ui\ui_manager.py                |    372 |    137 |   63.2% |
| moleditpy\src\moleditpy\ui\user_template_dialog.py      |    366 |    155 |   57.7% |
| moleditpy\src\moleditpy\ui\view_3d_logic.py             |    986 |    246 |   75.1% |
| moleditpy\src\moleditpy\ui\zoomable_view.py             |     72 |     37 |   48.6% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_2d_tab.py |     79 |     68 |   13.9% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_3d_tabs.py |    135 |    119 |   11.9% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_other_tab.py |     52 |     42 |   19.2% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_tab_base.py |     34 |     24 |   29.4% |
| moleditpy\src\moleditpy\utils\constants.py              |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\default_settings.py       |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\sip_isdeleted_safe.py     |     19 |      8 |   57.9% |
| moleditpy\src\moleditpy\utils\system_utils.py           |     34 |     18 |   47.1% |
| **TOTAL** | **13734** | **4476** | **67.41%** |

## Test Suite Status
- **Unit tests**: PASSED
- **Integration tests**: PASSED
- **GUI tests**: PASSED

[View Detailed HTML Report](coverage_html/index.html)