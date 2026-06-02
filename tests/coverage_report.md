# MoleditPy Coverage Report

- **Overall Project Coverage (Full)**: **80.17%**
- **Core Molecular Logic Coverage**: **81.90%**

> [!NOTE]
> **Core Molecular Logic Coverage** excludes UI boilerplate (dialogs, view managers, and interactor styles) to focus on scientific algorithm reliability.

### Core Logic Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      6 |      2 |   66.7% |
| moleditpy\src\moleditpy\main.py                         |     86 |     28 |   67.4% |
| moleditpy\src\moleditpy\core\mol_geometry.py            |    260 |     19 |   92.7% |
| moleditpy\src\moleditpy\core\molecular_data.py          |    307 |     38 |   87.6% |
| moleditpy\src\moleditpy\plugins\plugin_interface.py     |    128 |      0 |  100.0% |
| moleditpy\src\moleditpy\plugins\plugin_manager.py       |    370 |     15 |   95.9% |
| moleditpy\src\moleditpy\plugins\plugin_manager_window.py |    183 |      2 |   98.9% |
| moleditpy\src\moleditpy\ui\about_dialog.py              |     66 |     15 |   77.3% |
| moleditpy\src\moleditpy\ui\align_plane_dialog.py        |    161 |     17 |   89.4% |
| moleditpy\src\moleditpy\ui\alignment_dialog.py          |    165 |     11 |   93.3% |
| moleditpy\src\moleditpy\ui\angle_dialog.py              |    255 |     31 |   87.8% |
| moleditpy\src\moleditpy\ui\app_state.py                 |    423 |     99 |   76.6% |
| moleditpy\src\moleditpy\ui\atom_item.py                 |    311 |     46 |   85.2% |
| moleditpy\src\moleditpy\ui\atom_picking.py              |    170 |     46 |   72.9% |
| moleditpy\src\moleditpy\ui\base_picking_dialog.py       |     77 |      6 |   92.2% |
| moleditpy\src\moleditpy\ui\bond_item.py                 |    377 |     64 |   83.0% |
| moleditpy\src\moleditpy\ui\bond_length_dialog.py        |    229 |     20 |   91.3% |
| moleditpy\src\moleditpy\ui\calculation_worker.py        |    550 |     89 |   83.8% |
| moleditpy\src\moleditpy\ui\color_settings_dialog.py     |    199 |     32 |   83.9% |
| moleditpy\src\moleditpy\ui\compute_logic.py             |    399 |     96 |   75.9% |
| moleditpy\src\moleditpy\ui\dialog_3d_picking_mixin.py   |    139 |     21 |   84.9% |
| moleditpy\src\moleditpy\ui\dialog_logic.py              |    253 |     43 |   83.0% |
| moleditpy\src\moleditpy\ui\dihedral_dialog.py           |    259 |     35 |   86.5% |
| moleditpy\src\moleditpy\ui\edit_3d_logic.py             |    246 |     41 |   83.3% |
| moleditpy\src\moleditpy\ui\edit_actions_logic.py        |    858 |    217 |   74.7% |
| moleditpy\src\moleditpy\ui\export_logic.py              |    534 |    160 |   70.0% |
| moleditpy\src\moleditpy\ui\io_logic.py                  |    637 |    131 |   79.4% |
| moleditpy\src\moleditpy\ui\main_window_init.py          |   1255 |    175 |   86.1% |
| moleditpy\src\moleditpy\ui\molecular_scene_handler.py   |    958 |    179 |   81.3% |
| moleditpy\src\moleditpy\ui\molecule_scene.py            |    533 |    131 |   75.4% |
| moleditpy\src\moleditpy\ui\move_group_dialog.py         |    384 |    135 |   64.8% |
| moleditpy\src\moleditpy\ui\move_selected_atoms_dialog.py |    396 |     46 |   88.4% |
| moleditpy\src\moleditpy\ui\settings_dialog.py           |    107 |     11 |   89.7% |
| moleditpy\src\moleditpy\ui\string_importers.py          |    128 |     27 |   78.9% |
| moleditpy\src\moleditpy\ui\template_preview_item.py     |    105 |      5 |   95.2% |
| moleditpy\src\moleditpy\ui\template_preview_view.py     |     45 |      4 |   91.1% |
| moleditpy\src\moleditpy\ui\ui_manager.py                |    388 |    118 |   69.6% |
| moleditpy\src\moleditpy\ui\user_template_dialog.py      |    365 |     85 |   76.7% |
| moleditpy\src\moleditpy\ui\view_3d_logic.py             |   1017 |    265 |   73.9% |
| moleditpy\src\moleditpy\ui\zoomable_view.py             |     77 |      0 |  100.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_2d_tab.py |    105 |      0 |  100.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_3d_tabs.py |    135 |      0 |  100.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_other_tab.py |     66 |      0 |  100.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_tab_base.py |     56 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\constants.py              |     54 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\default_settings.py       |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\sip_isdeleted_safe.py     |     19 |      0 |  100.0% |
| **TOTAL** | **13842** | **2505** | **81.90%** |

### Full Application Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      6 |      2 |   66.7% |
| moleditpy\src\moleditpy\main.py                         |     86 |     28 |   67.4% |
| moleditpy\src\moleditpy\core\mol_geometry.py            |    260 |     19 |   92.7% |
| moleditpy\src\moleditpy\core\molecular_data.py          |    307 |     38 |   87.6% |
| moleditpy\src\moleditpy\plugins\plugin_interface.py     |    128 |      0 |  100.0% |
| moleditpy\src\moleditpy\plugins\plugin_manager.py       |    370 |     15 |   95.9% |
| moleditpy\src\moleditpy\plugins\plugin_manager_window.py |    183 |      2 |   98.9% |
| moleditpy\src\moleditpy\ui\__init__.py                  |      7 |      3 |   57.1% |
| moleditpy\src\moleditpy\ui\about_dialog.py              |     66 |     15 |   77.3% |
| moleditpy\src\moleditpy\ui\align_plane_dialog.py        |    163 |     18 |   89.0% |
| moleditpy\src\moleditpy\ui\alignment_dialog.py          |    167 |     12 |   92.8% |
| moleditpy\src\moleditpy\ui\analysis_window.py           |    116 |     25 |   78.4% |
| moleditpy\src\moleditpy\ui\angle_dialog.py              |    257 |     32 |   87.5% |
| moleditpy\src\moleditpy\ui\app_state.py                 |    423 |     99 |   76.6% |
| moleditpy\src\moleditpy\ui\atom_item.py                 |    311 |     46 |   85.2% |
| moleditpy\src\moleditpy\ui\atom_picking.py              |    170 |     46 |   72.9% |
| moleditpy\src\moleditpy\ui\base_picking_dialog.py       |     79 |      7 |   91.1% |
| moleditpy\src\moleditpy\ui\bond_item.py                 |    377 |     64 |   83.0% |
| moleditpy\src\moleditpy\ui\bond_length_dialog.py        |    231 |     21 |   90.9% |
| moleditpy\src\moleditpy\ui\calculation_worker.py        |    550 |     89 |   83.8% |
| moleditpy\src\moleditpy\ui\color_settings_dialog.py     |    199 |     32 |   83.9% |
| moleditpy\src\moleditpy\ui\compute_logic.py             |    399 |     96 |   75.9% |
| moleditpy\src\moleditpy\ui\constrained_optimization_dialog.py |    471 |    120 |   74.5% |
| moleditpy\src\moleditpy\ui\custom_interactor_style.py   |    473 |    337 |   28.8% |
| moleditpy\src\moleditpy\ui\custom_qt_interactor.py      |     44 |     35 |   20.5% |
| moleditpy\src\moleditpy\ui\dialog_3d_picking_mixin.py   |    143 |     24 |   83.2% |
| moleditpy\src\moleditpy\ui\dialog_logic.py              |    253 |     43 |   83.0% |
| moleditpy\src\moleditpy\ui\dihedral_dialog.py           |    261 |     36 |   86.2% |
| moleditpy\src\moleditpy\ui\edit_3d_logic.py             |    246 |     41 |   83.3% |
| moleditpy\src\moleditpy\ui\edit_actions_logic.py        |    858 |    217 |   74.7% |
| moleditpy\src\moleditpy\ui\export_logic.py              |    534 |    160 |   70.0% |
| moleditpy\src\moleditpy\ui\geometry_base_dialog.py      |     61 |      6 |   90.2% |
| moleditpy\src\moleditpy\ui\io_logic.py                  |    637 |    131 |   79.4% |
| moleditpy\src\moleditpy\ui\main_window.py               |     76 |     19 |   75.0% |
| moleditpy\src\moleditpy\ui\main_window_init.py          |   1255 |    175 |   86.1% |
| moleditpy\src\moleditpy\ui\mirror_dialog.py             |     71 |      7 |   90.1% |
| moleditpy\src\moleditpy\ui\molecular_scene_handler.py   |    958 |    179 |   81.3% |
| moleditpy\src\moleditpy\ui\molecule_scene.py            |    533 |    131 |   75.4% |
| moleditpy\src\moleditpy\ui\move_group_dialog.py         |    384 |    135 |   64.8% |
| moleditpy\src\moleditpy\ui\move_selected_atoms_dialog.py |    396 |     46 |   88.4% |
| moleditpy\src\moleditpy\ui\periodic_table_dialog.py     |     34 |      4 |   88.2% |
| moleditpy\src\moleditpy\ui\planarize_dialog.py          |    111 |     13 |   88.3% |
| moleditpy\src\moleditpy\ui\settings_dialog.py           |    107 |     11 |   89.7% |
| moleditpy\src\moleditpy\ui\string_importers.py          |    128 |     27 |   78.9% |
| moleditpy\src\moleditpy\ui\template_preview_item.py     |    105 |      5 |   95.2% |
| moleditpy\src\moleditpy\ui\template_preview_view.py     |     45 |      4 |   91.1% |
| moleditpy\src\moleditpy\ui\translation_dialog.py        |    243 |     11 |   95.5% |
| moleditpy\src\moleditpy\ui\ui_manager.py                |    388 |    118 |   69.6% |
| moleditpy\src\moleditpy\ui\user_template_dialog.py      |    365 |     85 |   76.7% |
| moleditpy\src\moleditpy\ui\view_3d_logic.py             |   1017 |    265 |   73.9% |
| moleditpy\src\moleditpy\ui\zoomable_view.py             |     77 |      0 |  100.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_2d_tab.py |    105 |      0 |  100.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_3d_tabs.py |    135 |      0 |  100.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_other_tab.py |     66 |      0 |  100.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_tab_base.py |     58 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\constants.py              |     54 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\default_settings.py       |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\sip_isdeleted_safe.py     |     19 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\system_utils.py           |     36 |      0 |  100.0% |
| **TOTAL** | **15603** | **3094** | **80.17%** |

## Test Suite Status
- **Unit tests**: PASSED
- **Integration tests**: PASSED
- **GUI tests**: PASSED

[View Detailed HTML Report](coverage_html/index.html)