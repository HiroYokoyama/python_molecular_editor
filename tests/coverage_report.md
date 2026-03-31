# MoleditPy Coverage Report

- **Overall Project Coverage (Full)**: **63.37%**
- **Core Molecular Logic Coverage**: **81.66%**

> [!NOTE]
> **Core Molecular Logic Coverage** excludes UI boilerplate (dialogs, view managers, and interactor styles) to focus on scientific algorithm reliability.

### Core Logic Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      6 |      2 |   66.7% |
| moleditpy\src\moleditpy\core\mol_geometry.py            |    261 |     19 |   92.7% |
| moleditpy\src\moleditpy\core\molecular_data.py          |    302 |     55 |   81.8% |
| moleditpy\src\moleditpy\plugins\plugin_interface.py     |     96 |     17 |   82.3% |
| moleditpy\src\moleditpy\plugins\plugin_manager.py       |    325 |    108 |   66.8% |
| moleditpy\src\moleditpy\ui\calculation_worker.py        |    546 |     91 |   83.3% |
| moleditpy\src\moleditpy\ui\string_importers.py          |    166 |     26 |   84.3% |
| moleditpy\src\moleditpy\utils\constants.py              |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\default_settings.py       |      1 |      0 |  100.0% |
| **TOTAL** | **1734** | **318** | **81.66%** |

### Full Application Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      6 |      2 |   66.7% |
| moleditpy\src\moleditpy\core\mol_geometry.py            |    261 |     19 |   92.7% |
| moleditpy\src\moleditpy\core\molecular_data.py          |    302 |     55 |   81.8% |
| moleditpy\src\moleditpy\plugins\plugin_interface.py     |     96 |     17 |   82.3% |
| moleditpy\src\moleditpy\plugins\plugin_manager.py       |    325 |    108 |   66.8% |
| moleditpy\src\moleditpy\ui\__init__.py                  |      7 |      3 |   57.1% |
| moleditpy\src\moleditpy\ui\about_dialog.py              |     64 |     15 |   76.6% |
| moleditpy\src\moleditpy\ui\align_plane_dialog.py        |    125 |     42 |   66.4% |
| moleditpy\src\moleditpy\ui\alignment_dialog.py          |    143 |     47 |   67.1% |
| moleditpy\src\moleditpy\ui\analysis_window.py           |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\ui\angle_dialog.py              |    229 |     94 |   59.0% |
| moleditpy\src\moleditpy\ui\app_state.py                 |    421 |    102 |   75.8% |
| moleditpy\src\moleditpy\ui\atom_item.py                 |    263 |     32 |   87.8% |
| moleditpy\src\moleditpy\ui\base_picking_dialog.py       |     63 |     19 |   69.8% |
| moleditpy\src\moleditpy\ui\bond_item.py                 |    304 |     71 |   76.6% |
| moleditpy\src\moleditpy\ui\bond_length_dialog.py        |    208 |     87 |   58.2% |
| moleditpy\src\moleditpy\ui\calculation_worker.py        |    546 |     91 |   83.3% |
| moleditpy\src\moleditpy\ui\color_settings_dialog.py     |    197 |     65 |   67.0% |
| moleditpy\src\moleditpy\ui\compute_logic.py             |    401 |     96 |   76.1% |
| moleditpy\src\moleditpy\ui\constrained_optimization_dialog.py |    406 |    307 |   24.4% |
| moleditpy\src\moleditpy\ui\custom_interactor_style.py   |    456 |    392 |   14.0% |
| moleditpy\src\moleditpy\ui\custom_qt_interactor.py      |     43 |     35 |   18.6% |
| moleditpy\src\moleditpy\ui\dialog_3d_picking_mixin.py   |     95 |     52 |   45.3% |
| moleditpy\src\moleditpy\ui\dialog_logic.py              |    232 |    140 |   39.7% |
| moleditpy\src\moleditpy\ui\dihedral_dialog.py           |    222 |     88 |   60.4% |
| moleditpy\src\moleditpy\ui\edit_3d_logic.py             |    236 |    128 |   45.8% |
| moleditpy\src\moleditpy\ui\edit_actions_logic.py        |    852 |    376 |   55.9% |
| moleditpy\src\moleditpy\ui\export_logic.py              |    511 |    150 |   70.6% |
| moleditpy\src\moleditpy\ui\geometry_base_dialog.py      |     58 |     41 |   29.3% |
| moleditpy\src\moleditpy\ui\io_logic.py                  |    621 |    201 |   67.6% |
| moleditpy\src\moleditpy\ui\main_window.py               |     62 |     19 |   69.4% |
| moleditpy\src\moleditpy\ui\main_window_init.py          |   1180 |    165 |   86.0% |
| moleditpy\src\moleditpy\ui\mirror_dialog.py             |     70 |      7 |   90.0% |
| moleditpy\src\moleditpy\ui\molecular_scene_handler.py   |    860 |    180 |   79.1% |
| moleditpy\src\moleditpy\ui\molecule_scene.py            |    506 |    187 |   63.0% |
| moleditpy\src\moleditpy\ui\move_group_dialog.py         |    351 |    208 |   40.7% |
| moleditpy\src\moleditpy\ui\periodic_table_dialog.py     |     33 |      7 |   78.8% |
| moleditpy\src\moleditpy\ui\planarize_dialog.py          |    105 |     28 |   73.3% |
| moleditpy\src\moleditpy\ui\settings_dialog.py           |    108 |     90 |   16.7% |
| moleditpy\src\moleditpy\ui\sip_isdeleted_safe.py        |     13 |      4 |   69.2% |
| moleditpy\src\moleditpy\ui\string_importers.py          |    166 |     26 |   84.3% |
| moleditpy\src\moleditpy\ui\template_preview_item.py     |    101 |     77 |   23.8% |
| moleditpy\src\moleditpy\ui\template_preview_view.py     |     42 |     14 |   66.7% |
| moleditpy\src\moleditpy\ui\translation_dialog.py        |    110 |     35 |   68.2% |
| moleditpy\src\moleditpy\ui\ui_manager.py                |    335 |     93 |   72.2% |
| moleditpy\src\moleditpy\ui\user_template_dialog.py      |    366 |    130 |   64.5% |
| moleditpy\src\moleditpy\ui\view_3d_logic.py             |    986 |    516 |   47.7% |
| moleditpy\src\moleditpy\ui\zoomable_view.py             |     72 |     39 |   45.8% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_2d_tab.py |     80 |     68 |   15.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_3d_tabs.py |    136 |    119 |   12.5% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_other_tab.py |     52 |     42 |   19.2% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_tab_base.py |     34 |     24 |   29.4% |
| moleditpy\src\moleditpy\utils\constants.py              |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\default_settings.py       |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\sip_isdeleted_safe.py     |     19 |      8 |   57.9% |
| moleditpy\src\moleditpy\utils\system_utils.py           |     34 |     18 |   47.1% |
| **TOTAL** | **13660** | **5004** | **63.37%** |

## Test Suite Status
- **Unit tests**: PASSED
- **Integration tests**: PASSED
- **GUI tests**: PASSED

[View Detailed HTML Report](coverage_html/index.html)