# MoleditPy Coverage Report

- **Overall Project Coverage (Full)**: **64.89%**
- **Core Molecular Logic Coverage**: **80.01%**

> [!NOTE]
> **Core Molecular Logic Coverage** excludes UI boilerplate (dialogs, view managers, and interactor styles) to focus on scientific algorithm reliability.

### Core Logic Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      6 |      2 |   66.7% |
| moleditpy\src\moleditpy\core\mol_geometry.py            |    240 |     15 |   93.8% |
| moleditpy\src\moleditpy\core\molecular_data.py          |    301 |     55 |   81.7% |
| moleditpy\src\moleditpy\plugins\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\plugins\plugin_manager.py       |    274 |     65 |   76.3% |
| moleditpy\src\moleditpy\ui\calculation_worker.py        |    535 |     86 |   83.9% |
| moleditpy\src\moleditpy\ui\compute_engine.py            |    482 |    119 |   75.3% |
| moleditpy\src\moleditpy\ui\molecular_parsers.py         |    355 |     88 |   75.2% |
| moleditpy\src\moleditpy\ui\project_io.py                |    233 |     79 |   66.1% |
| moleditpy\src\moleditpy\ui\string_importers.py          |    163 |     26 |   84.0% |
| moleditpy\src\moleditpy\utils\constants.py              |     31 |      0 |  100.0% |
| **TOTAL** | **2677** | **535** | **80.01%** |

### Full Application Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      6 |      2 |   66.7% |
| moleditpy\src\moleditpy\core\mol_geometry.py            |    240 |     15 |   93.8% |
| moleditpy\src\moleditpy\core\molecular_data.py          |    301 |     55 |   81.7% |
| moleditpy\src\moleditpy\plugins\plugin_interface.py     |     57 |      0 |  100.0% |
| moleditpy\src\moleditpy\plugins\plugin_manager.py       |    274 |     65 |   76.3% |
| moleditpy\src\moleditpy\ui\__init__.py                  |      7 |      3 |   57.1% |
| moleditpy\src\moleditpy\ui\about_dialog.py              |     64 |     15 |   76.6% |
| moleditpy\src\moleditpy\ui\align_plane_dialog.py        |    162 |     61 |   62.3% |
| moleditpy\src\moleditpy\ui\alignment_dialog.py          |    142 |     47 |   66.9% |
| moleditpy\src\moleditpy\ui\analysis_window.py           |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\ui\angle_dialog.py              |    297 |    144 |   51.5% |
| moleditpy\src\moleditpy\ui\app_state.py                 |    464 |    102 |   78.0% |
| moleditpy\src\moleditpy\ui\atom_item.py                 |    264 |     34 |   87.1% |
| moleditpy\src\moleditpy\ui\bond_item.py                 |    317 |     78 |   75.4% |
| moleditpy\src\moleditpy\ui\bond_length_dialog.py        |    269 |    136 |   49.4% |
| moleditpy\src\moleditpy\ui\calculation_worker.py        |    535 |     86 |   83.9% |
| moleditpy\src\moleditpy\ui\color_settings_dialog.py     |    181 |     53 |   70.7% |
| moleditpy\src\moleditpy\ui\compute_engine.py            |    482 |    119 |   75.3% |
| moleditpy\src\moleditpy\ui\compute_logic.py             |    344 |    147 |   57.3% |
| moleditpy\src\moleditpy\ui\constrained_optimization_dialog.py |    406 |    307 |   24.4% |
| moleditpy\src\moleditpy\ui\custom_interactor_style.py   |    465 |    401 |   13.8% |
| moleditpy\src\moleditpy\ui\custom_qt_interactor.py      |     43 |     35 |   18.6% |
| moleditpy\src\moleditpy\ui\dialog_3d_picking_mixin.py   |     92 |     50 |   45.7% |
| moleditpy\src\moleditpy\ui\dialog_logic.py              |    240 |    157 |   34.6% |
| moleditpy\src\moleditpy\ui\dihedral_dialog.py           |    300 |    139 |   53.7% |
| moleditpy\src\moleditpy\ui\edit_3d_logic.py             |    235 |    123 |   47.7% |
| moleditpy\src\moleditpy\ui\edit_actions_logic.py        |    775 |    359 |   53.7% |
| moleditpy\src\moleditpy\ui\export_logic.py              |    541 |    162 |   70.1% |
| moleditpy\src\moleditpy\ui\main_window.py               |    108 |     19 |   82.4% |
| moleditpy\src\moleditpy\ui\main_window_init.py          |   1158 |    151 |   87.0% |
| moleditpy\src\moleditpy\ui\mirror_dialog.py             |     70 |      7 |   90.0% |
| moleditpy\src\moleditpy\ui\molecular_parsers.py         |    355 |     88 |   75.2% |
| moleditpy\src\moleditpy\ui\molecular_scene_handler.py   |    857 |    170 |   80.2% |
| moleditpy\src\moleditpy\ui\molecule_scene.py            |    489 |    173 |   64.6% |
| moleditpy\src\moleditpy\ui\move_group_dialog.py         |    376 |    215 |   42.8% |
| moleditpy\src\moleditpy\ui\periodic_table_dialog.py     |     33 |      7 |   78.8% |
| moleditpy\src\moleditpy\ui\planarize_dialog.py          |    129 |     37 |   71.3% |
| moleditpy\src\moleditpy\ui\project_io.py                |    233 |     79 |   66.1% |
| moleditpy\src\moleditpy\ui\settings_dialog.py           |     98 |     38 |   61.2% |
| moleditpy\src\moleditpy\ui\sip_isdeleted_safe.py        |     13 |      4 |   69.2% |
| moleditpy\src\moleditpy\ui\string_importers.py          |    163 |     26 |   84.0% |
| moleditpy\src\moleditpy\ui\template_preview_item.py     |    101 |     77 |   23.8% |
| moleditpy\src\moleditpy\ui\template_preview_view.py     |     42 |     14 |   66.7% |
| moleditpy\src\moleditpy\ui\translation_dialog.py        |    197 |     87 |   55.8% |
| moleditpy\src\moleditpy\ui\ui_manager.py                |    305 |     84 |   72.5% |
| moleditpy\src\moleditpy\ui\user_template_dialog.py      |    370 |    131 |   64.6% |
| moleditpy\src\moleditpy\ui\view_3d_logic.py             |    930 |    480 |   48.4% |
| moleditpy\src\moleditpy\ui\view_loaders.py              |    155 |    104 |   32.9% |
| moleditpy\src\moleditpy\ui\zoomable_view.py             |     72 |     39 |   45.8% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_2d_tab.py |    100 |      9 |   91.0% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_3d_tabs.py |    174 |     26 |   85.1% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_other_tab.py |     42 |      1 |   97.6% |
| moleditpy\src\moleditpy\ui\settings_tabs\settings_tab_base.py |     11 |      3 |   72.7% |
| moleditpy\src\moleditpy\utils\constants.py              |     31 |      0 |  100.0% |
| moleditpy\src\moleditpy\utils\sip_isdeleted_safe.py     |     19 |      8 |   57.9% |
| moleditpy\src\moleditpy\utils\system_utils.py           |     34 |     18 |   47.1% |
| **TOTAL** | **14282** | **5015** | **64.89%** |

## Test Suite Status
- **Unit tests**: PASSED
- **Integration tests**: PASSED
- **GUI tests**: PASSED

[View Detailed HTML Report](coverage_html/index.html)