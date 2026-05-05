# MoleditPy Coverage Report

- **Overall Project Coverage (Full)**: **78.82%**
- **Core Molecular Logic Coverage**: **81.30%**

> [!NOTE]
> **Core Molecular Logic Coverage** excludes UI boilerplate (dialogs, view managers, and interactor styles) to focus on scientific algorithm reliability.

### Core Logic Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      6 |      2 |   66.7% |
| moleditpy\src\moleditpy\main.py                         |     80 |     22 |   72.5% |
| moleditpy\src\moleditpy\core\mol_geometry.py            |    260 |     19 |   92.7% |
| moleditpy\src\moleditpy\core\molecular_data.py          |    307 |     38 |   87.6% |
| moleditpy\src\moleditpy\plugins\plugin_interface.py     |    128 |      0 |  100.0% |
| moleditpy\src\moleditpy\plugins\plugin_manager.py       |    370 |     15 |   95.9% |
| moleditpy\src\moleditpy\plugins\plugin_manager_window.py |    183 |      2 |   98.9% |
| moleditpy\src\moleditpy\ui\about_dialog.py              |     66 |     15 |   77.3% |
| moleditpy\src\moleditpy\ui\align_plane_dialog.py        |    123 |     17 |   86.2% |
| moleditpy\src\moleditpy\ui\alignment_dialog.py          |    140 |     13 |   90.7% |
| moleditpy\src\moleditpy\ui\angle_dialog.py              |    233 |     20 |   91.4% |
| moleditpy\src\moleditpy\ui\app_state.py                 |    423 |     99 |   76.6% |
| moleditpy\src\moleditpy\ui\atom_item.py                 |    267 |     35 |   86.9% |
| moleditpy\src\moleditpy\ui\base_picking_dialog.py       |     70 |      5 |   92.9% |
| moleditpy\src\moleditpy\ui\bond_item.py                 |    301 |     71 |   76.4% |
| moleditpy\src\moleditpy\ui\bond_length_dialog.py        |    211 |     14 |   93.4% |
| moleditpy\src\moleditpy\ui\calculation_worker.py        |    549 |     92 |   83.2% |
| moleditpy\src\moleditpy\ui\color_settings_dialog.py     |    199 |     32 |   83.9% |
| moleditpy\src\moleditpy\ui\compute_logic.py             |    399 |     96 |   75.9% |
| moleditpy\src\moleditpy\ui\dialog_3d_picking_mixin.py   |    138 |     21 |   84.8% |
| moleditpy\src\moleditpy\ui\dialog_logic.py              |    228 |     31 |   86.4% |
| moleditpy\src\moleditpy\ui\dihedral_dialog.py           |    231 |     19 |   91.8% |
| moleditpy\src\moleditpy\ui\edit_3d_logic.py             |    243 |     44 |   81.9% |
| moleditpy\src\moleditpy\ui\edit_actions_logic.py        |    858 |    217 |   74.7% |
| moleditpy\src\moleditpy\ui\export_logic.py              |    512 |    150 |   70.7% |
| moleditpy\src\moleditpy\ui\io_logic.py                  |    637 |    131 |   79.4% |
| moleditpy\src\moleditpy\ui\main_window_init.py          |   1227 |    169 |   86.2% |
| moleditpy\src\moleditpy\ui\molecular_scene_handler.py   |    859 |    178 |   79.3% |
| moleditpy\src\moleditpy\ui\molecule_scene.py            |    502 |    145 |   71.1% |
| moleditpy\src\moleditpy\ui\move_group_dialog.py         |    399 |    159 |   60.2% |
| moleditpy\src\moleditpy\ui\settings_dialog.py           |    107 |     11 |   89.7% |
| moleditpy\src\moleditpy\ui\string_importers.py          |    128 |     27 |   78.9% |
| moleditpy\src\moleditpy\ui\template_preview_item.py     |    105 |      5 |   95.2% |
| moleditpy\src\moleditpy\ui\template_preview_view.py     |     45 |      4 |   91.1% |
| moleditpy\src\moleditpy\ui\ui_manager.py                |    382 |    116 |   69.6% |
| moleditpy\src\moleditpy\ui\user_template_dialog.py      |    365 |     85 |   76.7% |
| moleditpy\src\moleditpy\ui\view_3d_logic.py             |   1002 |    252 |   74.9% |
| moleditpy\src\moleditpy\ui\zoomable_view.py             |     77 |      0 |  100.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_2d_tab.py |     81 |      0 |  100.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_3d_tabs.py |    135 |      0 |  100.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_other_tab.py |     54 |      0 |  100.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_tab_base.py |     34 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\constants.py              |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\default_settings.py       |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\sip_isdeleted_safe.py     |     19 |      7 |   63.2% |
| **TOTAL** | **12715** | **2378** | **81.30%** |

### Full Application Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      6 |      2 |   66.7% |
| moleditpy\src\moleditpy\main.py                         |     80 |     22 |   72.5% |
| moleditpy\src\moleditpy\core\mol_geometry.py            |    260 |     19 |   92.7% |
| moleditpy\src\moleditpy\core\molecular_data.py          |    307 |     38 |   87.6% |
| moleditpy\src\moleditpy\plugins\plugin_interface.py     |    128 |      0 |  100.0% |
| moleditpy\src\moleditpy\plugins\plugin_manager.py       |    370 |     15 |   95.9% |
| moleditpy\src\moleditpy\plugins\plugin_manager_window.py |    183 |      2 |   98.9% |
| moleditpy\src\moleditpy\ui\__init__.py                  |      7 |      3 |   57.1% |
| moleditpy\src\moleditpy\ui\about_dialog.py              |     66 |     15 |   77.3% |
| moleditpy\src\moleditpy\ui\align_plane_dialog.py        |    125 |     18 |   85.6% |
| moleditpy\src\moleditpy\ui\alignment_dialog.py          |    142 |     14 |   90.1% |
| moleditpy\src\moleditpy\ui\analysis_window.py           |    116 |     25 |   78.4% |
| moleditpy\src\moleditpy\ui\angle_dialog.py              |    235 |     21 |   91.1% |
| moleditpy\src\moleditpy\ui\app_state.py                 |    425 |    100 |   76.5% |
| moleditpy\src\moleditpy\ui\atom_item.py                 |    267 |     35 |   86.9% |
| moleditpy\src\moleditpy\ui\base_picking_dialog.py       |     72 |      6 |   91.7% |
| moleditpy\src\moleditpy\ui\bond_item.py                 |    301 |     71 |   76.4% |
| moleditpy\src\moleditpy\ui\bond_length_dialog.py        |    213 |     15 |   93.0% |
| moleditpy\src\moleditpy\ui\calculation_worker.py        |    549 |     92 |   83.2% |
| moleditpy\src\moleditpy\ui\color_settings_dialog.py     |    199 |     32 |   83.9% |
| moleditpy\src\moleditpy\ui\compute_logic.py             |    399 |     96 |   75.9% |
| moleditpy\src\moleditpy\ui\constrained_optimization_dialog.py |    442 |    144 |   67.4% |
| moleditpy\src\moleditpy\ui\custom_interactor_style.py   |    460 |    395 |   14.1% |
| moleditpy\src\moleditpy\ui\custom_qt_interactor.py      |     44 |     35 |   20.5% |
| moleditpy\src\moleditpy\ui\dialog_3d_picking_mixin.py   |    142 |     24 |   83.1% |
| moleditpy\src\moleditpy\ui\dialog_logic.py              |    228 |     31 |   86.4% |
| moleditpy\src\moleditpy\ui\dihedral_dialog.py           |    233 |     20 |   91.4% |
| moleditpy\src\moleditpy\ui\edit_3d_logic.py             |    243 |     44 |   81.9% |
| moleditpy\src\moleditpy\ui\edit_actions_logic.py        |    858 |    217 |   74.7% |
| moleditpy\src\moleditpy\ui\export_logic.py              |    512 |    150 |   70.7% |
| moleditpy\src\moleditpy\ui\geometry_base_dialog.py      |     61 |      8 |   86.9% |
| moleditpy\src\moleditpy\ui\io_logic.py                  |    637 |    131 |   79.4% |
| moleditpy\src\moleditpy\ui\main_window.py               |     76 |     19 |   75.0% |
| moleditpy\src\moleditpy\ui\main_window_init.py          |   1227 |    169 |   86.2% |
| moleditpy\src\moleditpy\ui\mirror_dialog.py             |     71 |      7 |   90.1% |
| moleditpy\src\moleditpy\ui\molecular_scene_handler.py   |    859 |    178 |   79.3% |
| moleditpy\src\moleditpy\ui\molecule_scene.py            |    502 |    145 |   71.1% |
| moleditpy\src\moleditpy\ui\move_group_dialog.py         |    399 |    159 |   60.2% |
| moleditpy\src\moleditpy\ui\periodic_table_dialog.py     |     34 |      7 |   79.4% |
| moleditpy\src\moleditpy\ui\planarize_dialog.py          |    110 |     13 |   88.2% |
| moleditpy\src\moleditpy\ui\settings_dialog.py           |    107 |     11 |   89.7% |
| moleditpy\src\moleditpy\ui\string_importers.py          |    128 |     27 |   78.9% |
| moleditpy\src\moleditpy\ui\template_preview_item.py     |    105 |      5 |   95.2% |
| moleditpy\src\moleditpy\ui\template_preview_view.py     |     45 |      4 |   91.1% |
| moleditpy\src\moleditpy\ui\translation_dialog.py        |    243 |     11 |   95.5% |
| moleditpy\src\moleditpy\ui\ui_manager.py                |    382 |    116 |   69.6% |
| moleditpy\src\moleditpy\ui\user_template_dialog.py      |    365 |     85 |   76.7% |
| moleditpy\src\moleditpy\ui\view_3d_logic.py             |   1002 |    252 |   74.9% |
| moleditpy\src\moleditpy\ui\zoomable_view.py             |     77 |      0 |  100.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_2d_tab.py |     81 |      0 |  100.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_3d_tabs.py |    135 |      0 |  100.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_other_tab.py |     54 |      0 |  100.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_tab_base.py |     36 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\constants.py              |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\default_settings.py       |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\sip_isdeleted_safe.py     |     19 |      7 |   63.2% |
| moleditpy\src\moleditpy\utils\system_utils.py           |     36 |      2 |   94.4% |
| **TOTAL** | **14435** | **3057** | **78.82%** |

## Test Suite Status
- **Unit tests**: PASSED
- **Integration tests**: PASSED
- **GUI tests**: PASSED

[View Detailed HTML Report](coverage_html/index.html)