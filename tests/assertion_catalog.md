# Test Assertions Catalog

## tests/unit/test_app_logic.py

### test_worker_halt_logic
_Verify that CalculationWorker respects the halt flag (software flow)._

- assert 'Halted' in str(blocker.args[0])

### test_ez_preservation_logic
_Verify worker preserves explicit E/Z labels even if RDKit might lose them._

- assert bond.GetStereo() == Chem.BondStereo.STEREOZ

### test_molecular_data_fallback_serialization
_Verify MolecularData uses manual string construction if RDKit fails (fallback logic)._

- assert mol_block is not None
- assert 'MoleditPy' in mol_block
- assert 'X  ' in mol_block
- assert 'C  ' in mol_block

### test_coordinate_mapping_primary
_Verify primary RDKit conversion logic uses 'pos' attribute correctly._

- assert len(atom_lines) == 2
- assert abs(x_coord - expected_x) < 0.0001

### test_coordinate_mapping_fallback
_Verify fallback serialization (reverted logic) uses atom['item'].pos()._

- assert mol_block is not None
- assert 'MoleditPy' in mol_block
- assert len(atom_lines) == 2
- assert abs(x_coord - expected_x) < 0.0001

### test_worker_direct_mode
_Verify worker handles 'direct' mode (software branching)._

- assert mol.GetNumAtoms() == 5

## tests/unit/test_app_state.py

### test_get_current_state_captures_atoms
_get_current_state should capture all atoms with correct properties._

- assert len(atoms) == 2
- assert symbols == {'C', 'O'}

### test_get_current_state_captures_bonds
_get_current_state should capture bonds with correct order._

- assert len(state['bonds']) == 1
- assert bond_data['order'] == 2

### test_get_current_state_with_3d_mol
_State should include 3D molecule binary when current_mol is set._

- assert 'mol_3d' in state
- assert state['mol_3d'] is not None
- assert 'mol_3d_atom_ids' in state

### test_get_current_state_includes_version
_State should always include version string._

- assert 'version' in state

### test_get_current_state_captures_constraints
_3D constraints should be serialized in JSON-safe format._

- assert len(state['constraints_3d']) == 2
- assert isinstance(c, list)

### test_json_roundtrip_preserves_atoms
_Atoms should survive JSON serialization round-trip with correct properties._

- assert len(app.data.atoms) == 2
- assert symbols == {'C', 'N'}
- assert len(charged) == 1
- assert charged[0]['charge'] == 1

### test_json_roundtrip_preserves_bonds
_Bond order and stereo should survive JSON round-trip._

- assert len(app.data.bonds) == 1
- assert bond['order'] == 2

### test_json_roundtrip_preserves_radical
_Radical electrons should survive JSON round-trip._

- assert len(app.data.atoms) == 1
- assert next(iter(app.data.atoms.values()))['radical'] == 2

### test_get_current_state_empty
_Empty state should produce valid structure with no atoms/bonds._

- assert state['atoms'] == {}
- assert state['bonds'] == {}

## tests/unit/test_calculation_worker_direct.py

### test_calculation_worker_init
_No description provided._

- assert worker is not None
- assert hasattr(worker, 'status_update')
- assert hasattr(worker, 'finished')
- assert hasattr(worker, 'error')

### test_calculation_worker_halt_logic
_Test the internal _check_halted logic via run_calculation._

- assert any(('Halted' in str(val) for val in error_captor.emitted_values))

### test_calculation_worker_direct_mode
_Test 'direct' conversion mode which avoids RDKit 3D embedding._

- assert len(finish_captor.emitted_values) > 0
- assert res_mol.GetNumAtoms() > 2
- assert any((z > 0 for z in z_coords))

### test_calculation_worker_explicit_stereo_m_cfg
_Test parsing of M CFG labels in MOL block._

- assert len(finish_captor.emitted_values) > 0
- assert bond.GetStereo() == Chem.BondStereo.STEREOE

### test_calculation_worker_error_empty_input
_No description provided._

- assert any(('No atoms to convert' in str(val) for val in error_captor.emitted_values))

### test_calculation_worker_safe_helpers_halted
_Test that safe helpers don't emit finished if halted._

- assert len(finish_captor.emitted_values) == 0

### test_calculation_worker_rdkit_embedding_fail_fallback
_Test that embedding failure triggers fallback status or error messages._

- assert any(('failed' in m or 'error' in m for m in all_msgs))

## tests/unit/test_compute_logic.py

### test_on_calculation_error_stale
_Test on_calculation_error when the worker is stale (not in active set)._

- assert not compute.statusBar().showMessage.called

### test_on_calculation_error_basic
_Test on_calculation_error for an ACTIVE worker._

- assert compute.statusBar().showMessage.called
- assert 'Real Error' in compute.statusBar().showMessage.call_args[0][0]

### test_compute_set_optimization_method
_No description provided._

- assert compute.settings['optimization_method'] == 'GAFF_OBABEL'
- assert compute.statusBar().showMessage.called

### test_compute_halt_logic
_No description provided._

- assert 'test_id' in compute.halt_ids
- assert len(compute.active_worker_ids) == 0
- assert compute.statusBar().showMessage.called

### test_on_calculation_finished_basic
_No description provided._

- assert compute.current_mol == mol
- assert worker_id not in compute.active_worker_ids

### test_check_chemistry_problems_fallback_detects
_No description provided._

- assert c_item.has_problem is True
- assert compute.statusBar().showMessage.called

### test_trigger_conversion_empty
_No description provided._

- assert compute.statusBar().showMessage.called

### test_trigger_conversion_with_atoms
_No description provided._

- assert compute.statusBar().showMessage.called

### test_optimize_3d_structure_logic
_No description provided._

- assert compute.statusBar().showMessage.called

### test_on_calculation_finished_worker_id_mismatch
_No description provided._

- assert compute.current_mol is None

### test_trigger_conversion_chemistry_problems
_Test trigger_conversion when Chem.DetectChemistryProblems finds issues._

- assert any(('chemistry problem(s) found' in msg for msg in all_messages))

### test_trigger_conversion_sanitize_error
_Test trigger_conversion when Chem.SanitizeMol fails._

- assert any(('Error: Invalid chemical structure.' in msg for msg in all_messages))

### test_trigger_conversion_multiple_frags
_Test trigger_conversion with multiple fragments._

- assert any(('collision detection' in msg for msg in all_messages))

### test_on_calculation_finished_single_mol_legacy
_Test on_calculation_finished with a single mol (legacy result format)._

- assert compute.current_mol == mol

### test_on_calculation_error_legacy_payload
_Test on_calculation_error with a string (legacy error format)._

- assert any(('Fatal Error' in str(call[0][0]) for call in compute.statusBar().showMessage.call_args_list if call[0]))

## tests/unit/test_edit_3d_logic.py

### test_calculate_distance_logic
_No description provided._

- assert dist == pytest.approx(1.5)

### test_calculate_angle_logic
_No description provided._

- assert angle == pytest.approx(90.0)

### test_calculate_dihedral_logic
_No description provided._

- assert abs(dihedral) > 0

### test_handle_measurement_atom_selection
_No description provided._

- assert 10 in edit3d.selected_atoms_for_measurement
- assert len(edit3d.selected_atoms_for_measurement) == 1

### test_clear_3d_selection
_No description provided._

- assert edit3d.plotter.remove_actor.called

### test_toggle_measurement_mode
_No description provided._

- assert edit3d.measurement_mode is True
- assert edit3d.statusBar().showMessage.called
- assert edit3d.measurement_mode is False

### test_calculate_and_display_measurements_trigger
_No description provided._

- assert mock_display.called
- assert 'Distance' in mock_display.call_args[0][0][0]

## tests/unit/test_edit_actions.py

### test_rotate_molecule_2d_basic
_Test 2D rotation of the entire molecule._

- assert a1.setPos.called
- assert a2.setPos.called

### test_resolve_overlapping_groups_basic
_Test overlapping group resolution._

- assert a2.setPos.called

### test_update_implicit_hydrogens_main_logic
_Test calculation of implicit hydrogens for display._

- assert a_item.implicit_h_count == 4

### test_clipboard_copy_serialization
_Test copy selection MimeData generation._

- assert mock_clipboard.setMimeData.called

## tests/unit/test_export_logic.py

### test_create_multi_material_obj_advanced
_No description provided._

- assert os.path.exists(obj_path)
- assert os.path.exists(mtl_path)

### test_export_2d_png_basic_trigger
_No description provided._


### test_export_2d_svg_trigger
_No description provided._

- assert os.path.exists(save_path)

### test_export_stl_error_no_mol
_No description provided._


### test_export_obj_mtl_error_no_mol
_No description provided._


### test_export_stl_success_trigger
_No description provided._

- assert mesh.save.called

### test_export_obj_mtl_success_trigger
_No description provided._

- assert mock_create.called

### test_export_3d_png_logic
_No description provided._

- assert mock_plotter.screenshot.called

### test_export_color_stl_logic
_No description provided._

- assert exporter.statusBar().showMessage.called

### test_export_from_3d_view_with_colors_complex_splitting
_Test the complex logic of splitting a mesh by per-vertex colors._

- assert len(res) == 2
- assert colors[0] == (0, 0, 255)
- assert colors[1] == (255, 0, 0)

### test_export_2d_png_hides_items
_Test that export_2d_png hides non-atom items and restores them._

- assert other_item.hide.called

## tests/unit/test_geometry.py

### test_3d_bond_lengths
_Verify optimized 3D coordinates yield physical bond lengths._

- assert mol is not None
- assert cc_bond is not None
- assert 1.5 < dist < 1.6

### test_mirror_dialog_logic
_Verify MirrorDialog correctly manipulates coordinates and UI (software logic)._

- assert x_after == pytest.approx(-x_before, abs=0.001)
- assert main_window.draw_molecule_3d.called
- assert main_window.update_chiral_labels.called
- assert main_window.push_undo_state.called

### test_planarize_logic
_Verify planarize functionality using the actual PlanarizeDialog logic._

- assert s[-1] < 1e-10
- assert main_window.draw_molecule_3d.called
- assert main_window.update_chiral_labels.called
- assert main_window.push_undo_state.called

## tests/unit/test_hydrogen.py

### test_add_hydrogen_atoms_app_logic
_Verify add_hydrogen_atoms creates items in the scene based on RDKit logic._

- assert len(h_calls) == 4
- assert actions.scene.create_bond.call_count == 4

### test_remove_hydrogen_atoms_app_logic
_Verify remove_hydrogen_atoms finds and deletes H items using app logic._

- assert h_item in deleted_set
- assert actions.data.atoms[c_id]['item'] not in deleted_set

## tests/unit/test_io.py

### test_project_save_load_logic
_No description provided._

- assert os.path.exists(project_file)
- assert len(io_handler.data.atoms) == 1
- assert next(iter(io_handler.data.atoms.values()))['charge'] == 1

### test_xyz_export_logic
_No description provided._

- assert os.path.exists(xyz_file)
- assert xyz_mol.GetNumAtoms() == mol.GetNumAtoms()

## tests/unit/test_items_visual.py

### test_atom_item_visual_states
_Test AtomItem paint and boundingRect with different states._

- assert painter.drawEllipse.called or painter.drawText.called
- assert rect.width() > 0
- assert rect.height() > 0

### test_atom_item_hover
_Test AtomItem hover state changes._

- assert atom.hovered
- assert not atom.hovered

### test_atom_item_h_label_flip
_Test AtomItem flips H label position based on neighbor orientation._

- assert painter.drawText.called

### test_bond_item_render_complex_types
_Test BondItem rendering for triple bonds, wedge/dash, and E/Z labels._

- assert painter.drawLine.call_count == 3
- assert painter.drawPolygon.called
- assert painter.drawLine.call_count > 1
- assert painter.drawPath.called

### test_bond_item_ring_logic
_Test BondItem ring rendering logic by mocking RDKit mol integration._

- assert painter.drawLine.call_count >= 1

### test_atom_item_item_change_updates_bonds
_Test that moving an atom triggers bond position updates._

- assert bond.update_position.called

### test_atom_item_paint_transparent_bg
_Test AtomItem paint with transparent background uses CompositionMode_Clear._

- assert painter.setCompositionMode.called or painter.drawEllipse.called

### test_atom_item_paint_resilience_to_deleted_bond
_Test AtomItem paint doesn't crash when a C++ bond object is deleted._

- assert success

### test_atom_item_shape_collision
_Test AtomItem.shape() returns larger area than visual bounding for collision._

- assert isinstance(shape_path, QPainterPath)
- assert shape_rect.width() >= visual_rect.width() * 0.8
- assert shape_rect.height() >= visual_rect.height() * 0.8

### test_bond_item_shape_stroked
_Test bond shape is a stroked path (wider than line)._

- assert isinstance(path, QPainterPath)
- assert not path.isEmpty()
- assert rect.width() > 0
- assert rect.height() > 0

### test_bond_bounding_rect_geometry
_Test boundingRect includes bond line area properly._

- assert rect.width() >= 100
- assert rect.height() > 0

### test_bond_double_non_ring_parallel
_Test double bond parallel line rendering outside rings (non-ring case)._

- assert painter.drawLine.call_count == 2

### test_bond_hover_effects
_Test BondItem hover enter/leave effects by patching base classes to avoid type errors._

- assert bond.hovered
- assert scene.set_hovered_item.called
- assert not bond.hovered
- assert scene.set_hovered_item.call_count >= 2

### test_bond_stereo_e_rendering
_Test BondItem rendering for E-stereo label (stereo=4)._

- assert painter.drawPath.called

### test_bond_get_ez_label_rect
_Test get_ez_label_local_rect() calculates correct label geometry._

- assert label_rect is not None
- assert isinstance(label_rect, QRectF)
- assert label_rect.width() > 0

### test_bond_update_position_resilience
_Test BondItem.update_position resilience when atoms exist._

- assert True

## tests/unit/test_molecular_data.py

### test_add_atom_returns_incrementing_ids
_Verify atom IDs are sequential integers starting from 0._

- assert (id0, id1, id2) == (0, 1, 2)
- assert len(data.atoms) == 3

### test_add_atom_stores_properties
_Verify all atom properties (symbol, position, charge, radical) are stored._

- assert atom['symbol'] == 'N'
- assert atom['pos'].x() == pytest.approx(5.5)
- assert atom['pos'].y() == pytest.approx(-3.2)
- assert atom['charge'] == -1
- assert atom['radical'] == 1

### test_add_bond_creates_adjacency
_Verify bonds update the adjacency list bidirectionally._

- assert b in data.adjacency_list[a]
- assert a in data.adjacency_list[b]

### test_add_bond_normalizes_non_stereo_key
_Non-stereo bonds should have normalized keys (smaller ID first)._

- assert key == (a, b)
- assert status == 'created'

### test_add_bond_preserves_stereo_key_direction
_Stereo bonds should preserve the original key direction._

- assert key == (b, a)

### test_add_bond_update_existing
_Adding a bond with same atoms should update, not duplicate._

- assert status == 'updated'
- assert data.bonds[key]['order'] == 2

### test_remove_atom_cleans_adjacency_and_bonds
_Removing an atom should also remove its bonds and adjacency entries._

- assert b not in data.atoms
- assert b not in data.adjacency_list
- assert a not in data.adjacency_list.get(b, [])
- assert b not in key

### test_remove_bond_cleans_adjacency
_Removing a bond should update adjacency but leave atoms intact._

- assert len(data.bonds) == 0
- assert b not in data.adjacency_list[a]
- assert a not in data.adjacency_list[b]
- assert a in data.atoms and b in data.atoms

### test_remove_bond_reverse_lookup
_remove_bond should find bond regardless of key order._

- assert len(data.bonds) == 0

### test_to_rdkit_mol_ethanol
_Build ethanol in MolecularData and compare RDKit molecular formula._

- assert mol is not None
- assert rdMolDescriptors.CalcMolFormula(mol) == rdMolDescriptors.CalcMolFormula(ref_mol)

### test_to_rdkit_mol_benzene
_Build benzene (alternating single/double) and compare against RDKit reference._

- assert mol is not None
- assert mol.GetNumAtoms() == 6
- assert rdMolDescriptors.CalcMolFormula(mol) == rdMolDescriptors.CalcMolFormula(ref_mol)

### test_to_rdkit_mol_charged_atom
_Verify formal charges are correctly transferred to RDKit mol._

- assert mol is not None
- assert len(n_atoms) == 1
- assert n_atoms[0].GetFormalCharge() == 1

### test_to_rdkit_mol_radical
_Verify radical electrons are correctly transferred._

- assert mol is not None
- assert mol.GetAtomWithIdx(0).GetNumRadicalElectrons() == 1

### test_to_rdkit_mol_double_bond
_Verify double bond order is correctly transferred._

- assert mol is not None
- assert mol.GetBondBetweenAtoms(0, 1).GetBondTypeAsDouble() == 2.0

### test_to_rdkit_mol_triple_bond
_Verify triple bond: acetylene C#C._

- assert mol is not None
- assert mol.GetBondBetweenAtoms(0, 1).GetBondTypeAsDouble() == 3.0
- assert rdMolDescriptors.CalcMolFormula(mol) == rdMolDescriptors.CalcMolFormula(ref)

### test_to_rdkit_mol_coordinate_conversion
_Verify pixel-to-angstrom coordinate conversion is exact._

- assert pos.x == pytest.approx(px * ANGSTROM_PER_PIXEL, abs=1e-06)
- assert pos.y == pytest.approx(-py * ANGSTROM_PER_PIXEL, abs=1e-06)

### test_to_rdkit_mol_preserves_original_atom_id
_Verify _original_atom_id property is set on RDKit atoms._

- assert atom.HasProp('_original_atom_id')
- assert orig_id in [id0, id1]

### test_to_rdkit_mol_empty_returns_none
_Empty MolecularData should return None._

- assert data.to_rdkit_mol() is None

### test_to_mol_block_contains_all_atoms
_Verify MOL block text contains all expected atoms._

- assert mol_block is not None
- assert 'V2000' in mol_block or 'MoleditPy' in mol_block
- assert symbols == {'C', 'N', 'O'}

### test_molecular_weight_matches_rdkit
_Build acetic acid and compare MW against RDKit reference._

- assert Descriptors.HeavyAtomMolWt(mol) == pytest.approx(Descriptors.HeavyAtomMolWt(ref), abs=0.01)

## tests/unit/test_parsers.py

### test_fix_mol_block
_No description provided._

- assert 'V2000' in lines[3]
- assert len(lines[3]) >= 39

### test_load_mol_file_logic
_No description provided._

- assert len(parser.data.atoms) == 2
- assert any((a['symbol'] == 'C' for a in parser.data.atoms.values()))

### test_xyz_parsing_logic
_No description provided._

- assert mol is not None
- assert mol.GetNumAtoms() == 3
- assert any((a.GetSymbol() == 'O' for a in mol.GetAtoms()))

### test_load_xyz_file_with_estimation
_No description provided._

- assert mol is not None
- assert mol.GetNumAtoms() == 2
- assert mol.GetNumBonds() == 1

### test_save_as_mol_logic
_No description provided._

- assert os.path.exists(save_path)
- assert 'C  ' in content

## tests/unit/test_parsers_extended.py

### test_load_mol_file_fallback_to_sd_supplier
_No description provided._

- assert len(parser.data.atoms) == 1

### test_load_xyz_charge_prompt_ok
_No description provided._

- assert mol is not None

### test_load_xyz_skip_chemistry_via_button
_No description provided._

- assert mol is not None
- assert mol.HasProp('_xyz_skip_checks')

### test_load_xyz_unrecognized_symbol
_No description provided._


### test_load_xyz_unrecognized_symbol_auto_coerce
_No description provided._

- assert mol is not None
- assert mol.GetAtomWithIdx(0).GetSymbol() == 'C'

### test_save_as_xyz_logic
_No description provided._

- assert os.path.exists(save_path)

### test_load_mol_file_with_v2000_fix
_No description provided._

- assert len(parser.data.atoms) == 1

### test_load_xyz_file_with_manual_charge
_No description provided._

- assert mol is not None

### test_save_as_mol_logic
_No description provided._

- assert os.path.exists(save_path)

### test_estimate_bonds_from_distances
_No description provided._

- assert mol.GetNumBonds() == 1

### test_load_mol_file_cancel
_No description provided._

- assert not mock_dlg.called

### test_load_mol_file_dialog_cancel
_No description provided._

- assert parser.statusBar().showMessage.called or True

### test_load_mol_file_not_found
_No description provided._


### test_load_mol_file_invalid_format
_No description provided._

- assert any(('Invalid MOL file format' in str(call) for call in parser.statusBar().showMessage.call_args_list))

### test_load_xyz_file_invalid_atom_count
_No description provided._


### test_load_xyz_file_zero_atoms
_No description provided._


### test_load_xyz_file_too_few_lines
_No description provided._


### test_save_as_mol_cancel
_No description provided._


### test_save_as_xyz_cancel
_No description provided._


### test_fix_mol_counts_line
_No description provided._

- assert parser.fix_mol_counts_line(line1) == line1
- assert 'V2000' in parser.fix_mol_counts_line(line2)

### test_fix_mol_block
_No description provided._

- assert 'V2000' in fixed.splitlines()[3]

### test_load_xyz_always_ask_charge
_No description provided._

- assert mol.GetIntProp('_xyz_charge') == 1

### test_save_as_xyz_charge_mult
_No description provided._

- assert 'chrg = 1' in lines[1]

### test_load_mol_file_malformed_counts
_Test fix_mol_counts_line via load_mol_file._

- assert len(parser.data.atoms) == 1

### test_estimate_bonds_radius_fallback
_Test bond estimation with an unknown element radius in our dictionary._

- assert mol.GetNumBonds() == 1

### test_load_xyz_file_prompt_skip_logic
_Test the 'Skip chemistry' button branch in load_xyz_file._


### test_save_as_mol_no_current_path
_Test save_as_mol when no file is currently open (untitled)._

- assert os.path.exists(save_path)

### test_save_as_xyz_charge_exception
_Test save_as_xyz fallback when descriptors fail._

- assert os.path.exists(save_path)

## tests/unit/test_plugin_manager.py

### test_plugin_info_extracts_all_fields
_Verify AST parser extracts PLUGIN_NAME, VERSION, AUTHOR, DESCRIPTION._

- assert info['name'] == 'My Cool Plugin'
- assert info['version'] == '1.2.3'
- assert info['author'] == 'Test Author'
- assert info['description'] == 'Does amazing things'
- assert info['category'] == 'Analysis'

### test_plugin_info_version_tuple
_Version tuples like (1, 0, 0) should be joined as '1.0.0'._

- assert info['version'] == '2.3.1'

### test_plugin_info_fallback_dunder
___version__ and __author__ should be used when PLUGIN_ variants are missing._

- assert info['version'] == '0.5.0'
- assert info['author'] == 'Dunder Author'

### test_plugin_info_docstring_fallback
_Module docstring should be used as description if PLUGIN_DESCRIPTION is missing._

- assert info['description'] == 'This is the plugin description from docstring.'

### test_plugin_info_missing_file
_Non-existent file should return defaults, not crash._

- assert info['name'] == 'nonexistent.py'
- assert info['version'] == 'Unknown'

### test_plugin_info_syntax_error
_File with syntax errors should return defaults, not crash._

- assert info['version'] == 'Unknown'

### test_plugin_info_empty_file
_Empty file should return defaults._

- assert info['name'] == 'empty.py'

### test_register_menu_action
_register_menu_action should store action metadata._

- assert len(pm.menu_actions) == 1
- assert pm.menu_actions[0]['plugin'] == 'TestPlugin'

### test_register_toolbar_action
_No description provided._

- assert len(pm.toolbar_actions) == 1

### test_register_export_action
_No description provided._

- assert len(pm.export_actions) == 1
- assert pm.export_actions[0]['label'] == 'Export as PDF'

### test_register_optimization_method
_No description provided._

- assert 'MMFF94' in pm.optimization_methods

### test_register_file_opener
_No description provided._

- assert '.cif' in pm.file_openers

### test_register_file_opener_priority
_Higher priority opener should replace lower priority one in sorted order._

- assert pm.file_openers['.xyz'][0]['callback'] == cb_high

### test_register_analysis_tool
_No description provided._

- assert len(pm.analysis_tools) == 1

### test_register_save_handler
_No description provided._

- assert 'TestPlugin' in pm.save_handlers

### test_register_load_handler
_No description provided._

- assert 'TestPlugin' in pm.load_handlers

### test_register_3d_style
_No description provided._

- assert 'Wireframe' in pm.custom_3d_styles

### test_register_document_reset_handler
_No description provided._

- assert len(pm.document_reset_handlers) == 1

### test_invoke_document_reset_handlers
_All registered reset handlers should be called._

- assert called == ['P1', 'P2']

### test_discover_plugins_empty_dir
_Empty plugin directory should return no plugins._

- assert plugins == []

### test_discover_plugins_single_file
_Single .py file in plugin dir should be discovered._

- assert len(plugins) >= 1
- assert 'Hello' in names

### test_discover_plugins_package
_Folder with __init__.py should be discovered as package plugin._

- assert len(plugins) >= 1
- assert 'Package Plugin' in names

### test_discover_plugins_ignores_dunder_files
_Files starting with __ should be ignored, and root __init__.py shouldn't count as a plugin._

- assert len(plugins) == 1
- assert '__trash' not in names

### test_ensure_plugin_dir_creates_directory
_ensure_plugin_dir should create the directory if it doesn't exist._

- assert os.path.isdir(new_dir)

## tests/unit/test_project_io_extended.py

### test_save_project_no_data
_No description provided._


### test_save_project_overwrite_json
_No description provided._

- assert os.path.exists(project_file)
- assert data['format'] == 'PME Project'

### test_save_project_overwrite_raw
_No description provided._

- assert os.path.exists(raw_file)
- assert data['atoms'] == 'mock'

### test_save_project_redirect_to_save_as
_No description provided._

- assert mock_save_as.called

### test_load_raw_data_success
_No description provided._

- assert mock_set_state.called
- assert io.current_file_path == raw_file

### test_load_json_data_invalid_format
_No description provided._

- assert mock_warn.called

### test_open_project_file_dispatch
_No description provided._

- assert mock_json.called
- assert mock_raw.called

### test_save_as_json_trigger
_No description provided._

- assert os.path.exists(save_path)

### test_load_raw_data_error_paths
_No description provided._


## tests/unit/test_properties.py

### test_analysis_window_regular_mol
_Verify AnalysisWindow calculates and displays properties for regular molecules._

- assert formula_found
- assert 'C2H6O' in value

### test_analysis_window_xyz_derived
_Verify AnalysisWindow uses manual logic for XYZ-derived structures._

- assert 'C2HO' in formula_val
- assert not smiles_present

## tests/unit/test_scene_extended.py

### test_scene_keypress_modes
_No description provided._

- assert scene.window.activate_select_mode.called
- assert scene.mode == expected_mode

### test_scene_keypress_special_symbols
_No description provided._

- assert scene.mode == 'atom_Cl'
- assert scene.mode == 'atom_Br'
- assert scene.mode == 'atom_Si'

### test_scene_keypress_delete
_No description provided._


### test_scene_maintenance_methods
_No description provided._

- assert atom_item.has_problem is False

### test_scene_queries
_No description provided._

- assert scene.find_atom_near(QPointF(50, 0)) == a2

### test_scene_update_connected_bonds
_No description provided._

- assert a.bonds[0].update_position.called

### test_scene_leave_event
_No description provided._

- assert scene.template_preview.hide.called

### test_scene_update_bond_stereo
_No description provided._

- assert scene.data.bonds[10, 20]['stereo'] == 1
- assert bond.set_stereo.called

### test_scene_mouse_drag_create_bond_existing_atoms
_No description provided._

- assert (aid1, aid2) in scene.data.bonds or (aid2, aid1) in scene.data.bonds

### test_scene_mouse_click_create_single_atom
_No description provided._

- assert any((a['symbol'] == 'O' for a in scene.data.atoms.values()))

### test_scene_right_click_delete
_No description provided._


### test_scene_bond_complex_interactions
_No description provided._

- assert scene.data.bonds[bond_key]['stereo'] == 4
- assert scene.data.bonds[bond_key]['stereo'] == 3
- assert current_bond_data['order'] == 1
- assert current_bond_data['stereo'] == 1

### test_scene_atom_keyboard_properties
_No description provided._

- assert atom_item.charge == 1
- assert atom_item.charge == 0
- assert atom_item.radical == 1

### test_scene_template_interaction
_No description provided._


### test_scene_user_template_logic
_No description provided._

- assert any((a['symbol'] == 'N' and a['charge'] == 1 for a in scene.data.atoms.values()))

## tests/unit/test_scene_interactions.py

### test_scene_toggle_radical
_No description provided._

- assert atom_item.radical == 1
- assert atom_item.radical == 2
- assert atom_item.radical == 0

### test_scene_toggle_charge
_No description provided._

- assert atom_item.charge == 1
- assert atom_item.charge == 0

### test_add_benzene_fragment
_No description provided._

- assert len(scene.data.atoms) == 6
- assert len(scene.data.bonds) == 6
- assert len([b for b in scene.data.bonds.values() if b['order'] == 2]) == 3

### test_benzene_fusion_rotation
_No description provided._

- assert len(scene.data.atoms) == 6
- assert eb.order == 2

### test_delete_selected_items
_No description provided._


### test_double_click_select_component
_No description provided._

- assert scene.data.atoms[id1]['item'].setSelected.called
- assert scene.data.atoms[id2]['item'].setSelected.called
- assert scene.data.atoms[id3]['item'].setSelected.called
- assert not atom_iso.setSelected.called

### test_scene_key_event_dispatch
_Test key events like '4' (template), '.' (radical), +/- (charge)._

- assert atom_item.charge == 1
- assert scene.data.atoms[aid]['charge'] == 1

### test_scene_update_template_preview_logic
_Test update_template_preview with different targets._

- assert len(points) == 6
- assert 'items' in scene.template_context
- assert scene.template_context['items'][0] == atom_item

### test_scene_drag_create_bond_sequence
_Test the full mouse press -> move -> release sequence for creating a bond._

- assert bond is not None
- assert bond.order == 1
- assert getattr(scene, 'start_atom', None) == a1

## tests/unit/test_stereochemistry.py

### test_wedge_dash_mapping
_Verify that Wedge/Dash in MolecularData mapped to RDKit BondDir._

- assert mol is not None
- assert bond.GetBondDir() == Chem.BondDir.BEGINWEDGE
- assert bond.GetBondDir() == Chem.BondDir.ENDWEDGE
- assert bond.GetBondDir() == Chem.BondDir.BEGINDASH
- assert bond.GetBondDir() == Chem.BondDir.ENDDASH

### test_ez_stereo_persistence
_Verify E/Z double bond configurations are maintained._

- assert bond.GetStereo() == Chem.BondStereo.STEREOE
- assert bond_z.GetStereo() == Chem.BondStereo.STEREOZ

### test_chiral_r_s_consistency
_Verify that R/S markers match RDKit descriptors for a known chiral center._

- assert found_wedge

### test_stereo_confirmation
_Verify that mirroring a molecule actually inverts its stereochemistry._

- assert mol.GetAtomWithIdx(1).GetProp('_CIPCode') == 'S'
- assert mol.GetAtomWithIdx(1).GetProp('_CIPCode') == 'R'

### test_stereo_loss_on_planarize
_Verify that planarizing a chiral center removes its chirality._

- assert mol.GetAtomWithIdx(1).HasProp('_CIPCode')
- assert not mol.GetAtomWithIdx(1).HasProp('_CIPCode')

## tests/unit/test_string_importers.py

### test_smiles_ethanol_atom_count
_Ethanol SMILES should produce 3 heavy atoms matching RDKit reference._

- assert len(data.atoms) == ref_mol.GetNumAtoms()

### test_smiles_ethanol_bond_count
_Ethanol SMILES should produce correct bond count._

- assert len(data.bonds) == ref_mol.GetNumBonds()

### test_smiles_benzene
_Benzene SMILES should produce 6C atoms and 6 bonds._

- assert len(data.atoms) == 6
- assert len(data.bonds) == 6
- assert all((a['symbol'] == 'C' for a in data.atoms.values()))

### test_smiles_preserves_formal_charge
_Ammonium NH4+ should carry charge +1 on nitrogen._

- assert len(n_atoms) == 1
- assert n_atoms[0]['charge'] == 1

### test_smiles_aspirin_formula
_Aspirin SMILES should produce correct molecular formula vs RDKit._

- assert rdMolDescriptors.CalcMolFormula(mol) == rdMolDescriptors.CalcMolFormula(ref)

### test_smiles_stereocenter_wedge
_Chiral SMILES should produce at least one wedge or dash bond._

- assert len(stereo_bonds) >= 1

### test_smiles_invalid_shows_error
_Invalid SMILES should trigger status bar error, not crash._

- assert 'Invalid' in last_msg or 'Error' in last_msg

### test_smiles_empty_shows_error
_Empty SMILES should trigger error message._


### test_inchi_ethanol
_InChI import should match RDKit reference atom/bond count._

- assert len(data.atoms) == ref.GetNumAtoms()
- assert len(data.bonds) == ref.GetNumBonds()

### test_inchi_caffeine_formula
_Caffeine InChI should produce correct molecular formula._

- assert rdMolDescriptors.CalcMolFormula(mol) == rdMolDescriptors.CalcMolFormula(ref)

### test_inchi_invalid_shows_error
_Invalid InChI should show error, not crash._

- assert 'Invalid' in last_msg or 'Error' in last_msg

## tests/integration/test_calculation_worker.py

### test_calculation_worker_bond_length_validation
_Integration test: Run a real 3D optimization for Ethane._

- assert mol_3d is not None
- assert len(cc_bonds) >= 1
- assert dist_measured == pytest.approx(dist_ref, rel=0.01)

### test_calculation_worker_angle_validation
_Integration test: Run 3D optimization for Propane._

- assert central != -1
- assert angle_measured == pytest.approx(angle_ref, abs=5.0)

### test_calculation_worker_dihedral_validation
_Integration test: Run 3D optimization for n-Butane._

- assert is_staggered

### test_calculation_worker_conversion_no_optimize
_Integration test: 2D to 3D without optimization._

- assert 'F' in symbols
- assert 'Cl' in symbols
- assert dist > 0.1

### test_calculation_worker_direct_mode
_Integration test: Direct conversion mode._

- assert p1.z == pytest.approx(0.0)
- assert p2.z == pytest.approx(0.0)
- assert p1.x == pytest.approx(10.0 * scale)
- assert p1.y == pytest.approx(-20.0 * scale)
- assert len(h_atoms) > 0
- assert abs(hp.z) > 0.05

### test_calculation_worker_halt_logic
_Integration test: Verify halt mechanism._

- assert err_id == 123
- assert 'Halted' in err_msg

### test_calculation_worker_global_halt
_Integration test: Verify global halt._

- assert err_id is None
- assert 'Halted' in err_msg

### test_calculation_worker_invalid_input
_Test error handling for empty input._

- assert 'No atoms to convert' in err_msg

### test_calculation_worker_isolation
_Ensure worker_id correctly isolates halt signals._

- assert res_id == 789

### test_calculation_worker_direct_mode_stereo
_Integration test: Direct mode with wedge/dash bonds._

- assert p1.z == pytest.approx(0.0)
- assert p2.z == pytest.approx(0.0)
- assert p3.z == pytest.approx(1.5)

### test_calculation_worker_constraint_embedding_fallback
_Test the fallback to constraint-based embedding when initial embedding fails._

- assert mock_embed.call_count >= 2

### test_calculation_worker_uff_fallback
_Test fallback to UFF when MMFF optimization fails._

- assert mock_mmff.called
- assert mock_uff.called

### test_calculation_worker_mmff_variants
_Test switching between MMFF94 and MMFF94s._

- assert found

### test_calculation_worker_obabel_fallback_mocked
_Test the Open Babel fallback path by mocking availability and pybel._

- assert mock_pybel.readstring.called
- assert mock_ob_mol.make3D.called

### test_calculation_worker_complex_direct_h_placement
_Test direct mode with 4 hydrogens on a carbon to hit the rotation/offset logic._

- assert mol_3d.GetNumAtoms() >= 5

## tests/integration/test_chiral_labels_integration.py

### test_chiral_labels_toggle_3d
_Test that toggling 'Show Chiral Labels' correctly displays/hides labels in 3D._

- assert window.show_chiral_labels is False
- assert not labels_drawn()
- assert window.show_chiral_labels is True
- assert len(chiral_centers) == 1
- assert initial_label in ['R', 'S']
- assert initial_label in labels

### test_chiral_labels_mirror_inversion_3d
_Test that mirror transformation inverts the chiral label in 3D._

- assert initial_label in ['R', 'S']
- assert len(chiral_centers) == 1
- assert new_label != initial_label
- assert new_label in ['R', 'S']
- assert new_label in labels
- assert initial_label not in labels

## tests/gui/test_main_app.py

### test_molecular_data_add_atom
_MolecularData: 原子の追加テスト_

- assert atom_id == 0
- assert 0 in data.atoms
- assert data.atoms[0]['symbol'] == 'C'
- assert data.atoms[0]['pos'] == QPointF(0, 0)
- assert data._next_atom_id == 1
- assert 0 in data.adjacency_list

### test_molecular_data_add_bond
_MolecularData: 結合の追加テスト_

- assert key == (0, 1)
- assert status == 'created'
- assert (0, 1) in data.bonds
- assert data.bonds[0, 1]['order'] == 2
- assert data.bonds[0, 1]['stereo'] == 0
- assert data.adjacency_list[0] == [1]
- assert data.adjacency_list[1] == [0]

### test_molecular_data_add_stereo_bond
_MolecularData: 立体結合 (Wedge/Dash) がソートされずに保存されるかテスト_

- assert key == (1, 0)
- assert status == 'created'
- assert (1, 0) in data.bonds
- assert (0, 1) not in data.bonds
- assert data.bonds[1, 0]['stereo'] == 1

### test_molecular_data_remove_atom
_MolecularData: 原子削除と関連する結合の削除テスト_

- assert len(data.atoms) == 3
- assert len(data.bonds) == 2
- assert len(data.atoms) == 2
- assert 0 not in data.atoms
- assert 1 in data.atoms
- assert 2 in data.atoms
- assert len(data.bonds) == 0
- assert 0 not in data.adjacency_list
- assert 1 in data.adjacency_list
- assert data.adjacency_list[1] == []

### test_molecular_data_remove_bond
_MolecularData: 結合削除のテスト_

- assert len(data.bonds) == 1
- assert data.adjacency_list[0] == [1]
- assert len(data.bonds) == 0
- assert data.adjacency_list[0] == []
- assert data.adjacency_list[1] == []

### test_to_rdkit_mol_stereo
_MolecularData: 立体結合 (Wedge/Dash) のRDKit変換テスト_

- assert mol is not None
- assert wedge_bond.GetBondDir() == Chem.BondDir.BEGINWEDGE
- assert dash_bond.GetBondDir() == Chem.BondDir.BEGINDASH

### test_to_rdkit_mol_ez_stereo
_MolecularData: E/Z立体結合のRDKit変換テスト_

- assert mol is not None
- assert double_bond.GetBondType() == Chem.BondType.DOUBLE
- assert double_bond.GetStereo() == Chem.BondStereo.STEREOZ

### test_app_launch
_MainWindow: アプリケーションが正常に起動することを確認_

- assert window is not None
- assert window.isVisible()
- assert 'MoleditPy Ver.' in window.windowTitle()

### test_mode_change_atom
_ツールバー: 原子ボタンでモードが変更されることを確認_

- assert scene.mode == 'atom_C'
- assert n_button is not None
- assert scene.mode == 'atom_N'
- assert scene.current_atom_symbol == 'N'
- assert window.statusBar().currentMessage() == 'Mode: Draw Atom (N)'

### test_mode_change_bond
_ツールバー: 結合ボタンでモードが変更されることを確認_

- assert db_button is not None
- assert scene.mode == 'bond_2_0'
- assert scene.bond_order == 2
- assert scene.bond_stereo == 0
- assert window.statusBar().currentMessage() == 'Mode: Draw Bond (Order: 2)'

### test_draw_atom_on_click
_MoleculeScene: クリックで原子を描画するテスト_

- assert len(window.data.atoms) == 0
- assert len(window.data.atoms) == 1
- assert window.data.atoms[atom_id]['symbol'] == 'N'
- assert window.data.atoms[atom_id]['item'].pos() == click_pos

### test_draw_bond_on_drag
_MoleculeScene: ドラッグで結合を描画するテスト_

- assert len(window.data.atoms) == 0
- assert len(window.data.bonds) == 0
- assert len(window.data.atoms) == 2
- assert len(window.data.bonds) == 1
- assert (id1, id2) in window.data.bonds
- assert window.data.bonds[id1, id2]['order'] == 1

### test_draw_bond_to_existing_atom
_MoleculeScene: 既存の原子へドラッグして結合するテスト_

- assert len(window.data.atoms) == 2
- assert len(window.data.bonds) == 0
- assert len(window.data.atoms) == 2
- assert len(window.data.bonds) == 1
- assert (0, 1) in window.data.bonds

### test_change_atom_symbol_on_click
_MoleculeScene: 既存原子のクリックで元素を変更するテスト_

- assert window.data.atoms[0]['symbol'] == 'C'
- assert window.data.atoms[0]['symbol'] == 'O'
- assert len(window.data.atoms) == 1

### test_change_bond_order_on_click
_MoleculeScene: 既存結合のクリックで次数を変更するテスト_

- assert window.data.bonds[0, 1]['order'] == 1
- assert window.data.bonds[0, 1]['order'] == 2

### test_delete_atom_on_right_click
_MoleculeScene: 右クリックで原子を削除するテスト_

- assert len(window.data.atoms) == 1
- assert len(window.data.atoms) == 0

### test_charge_mode_click
_MoleculeScene: 電荷モードでのクリックテスト_

- assert window.data.atoms[0]['charge'] == 0
- assert window.data.atoms[0]['charge'] == 1
- assert window.data.atoms[0]['charge'] == -1

### test_2d_to_3d_conversion
_2D->3D変換: 変換ボタンのテスト_

- assert len(window.data.atoms) == 2
- assert len(window.data.bonds) == 1
- assert convert_button.isEnabled()
- assert window.current_mol is not None
- assert window.current_mol is not None
- assert window.optimize_3d_button.isEnabled()
- assert window.export_button.isEnabled()
- assert window.analysis_action.isEnabled()

### test_optimize_3d
_3D最適化: 3D最適化ボタンのテスト_

- assert window.current_mol is not None
- assert window.optimize_3d_button.isEnabled()
- assert 'optimization successful' in window.statusBar().currentMessage()

### test_change_3d_style
_3Dスタイル変更: スタイルメニューのテスト_

- assert window.current_3d_style == 'ball_and_stick'
- assert style_button is not None
- assert cpk_action is not None
- assert window.current_3d_style == 'cpk'

### test_undo_redo
_Undo/Redo: 操作のテスト_

- assert len(window.data.atoms) == 0
- assert isinstance(window.undo_stack, list)
- assert window.undo_action.isEnabled() in (True, False)
- assert window.redo_action.isEnabled() in (True, False)
- assert len(window.data.atoms) == 1
- assert window.undo_action.isEnabled() in (True, False)
- assert window.undo_action.isEnabled() in (True, False)
- assert window.redo_action.isEnabled() in (True, False)
- assert len(window.data.atoms) in (0, 1)
- assert window.undo_action.isEnabled() in (True, False)
- assert window.redo_action.isEnabled() in (True, False)
- assert len(window.data.atoms) == 1
- assert len(window.data.atoms) == 1
- assert window.undo_action.isEnabled() in (True, False)
- assert window.redo_action.isEnabled() in (True, False)

### test_clear_all
_Clear All: 全消去のテスト_

- assert len(window.data.atoms) == 1
- assert window.has_unsaved_changes == True
- assert len(window.data.atoms) == 0
- assert len(window.data.bonds) == 0
- assert window.current_mol is None
- assert window.has_unsaved_changes == False
- assert len(window.undo_stack) == 1

### test_copy_paste
_編集: コピー＆ペーストのテスト_

- assert len(window.data.atoms) == 2
- assert len(window.data.bonds) == 1
- assert len(scene.selectedItems()) == 3
- assert mime_data.hasFormat(moleditpy.CLIPBOARD_MIME_TYPE)
- assert len(window.data.atoms) == 4
- assert len(window.data.bonds) == 2
- assert window.data.atoms[2]['pos'].x() > 50
- assert window.data.atoms[2]['pos'].y() > 0
- assert window.data.atoms[3]['pos'].x() > 50
- assert window.data.atoms[3]['pos'].y() > 0

### test_file_import_smiles
_ファイル: SMILESインポートのテスト_

- assert len(window.data.atoms) == 3
- assert len(window.data.bonds) == 2
- assert symbols.count('C') == 2
- assert symbols.count('O') == 1
- assert 'Successfully loaded from SMILES' in window.statusBar().currentMessage()

### test_key_press_change_atom
_キーボードショートカット: 'O'キーで原子を変更_

- assert window.data.atoms[0]['symbol'] == 'C'
- assert window.data.atoms[0]['symbol'] == 'O'

### test_key_press_change_bond
_キーボードショートカット: '2'キーで結合次数を変更_

- assert window.data.bonds[0, 1]['order'] == 1
- assert window.data.bonds[0, 1]['order'] == 2

### test_radical_mode_toggle
_MoleculeScene: ラジカルモードでのクリックテスト_

- assert window.data.atoms[0]['radical'] == 0
- assert window.data.atoms[0]['radical'] == 1
- assert window.data.atoms[0]['radical'] == 2
- assert window.data.atoms[0]['radical'] == 0

### test_delete_key_selection
_MoleculeScene: Deleteキーで選択項目を削除_

- assert len(window.data.atoms) == 1
- assert len(scene.selectedItems()) == 1
- assert len(window.data.atoms) == 0

### test_draw_benzene_template
_MoleculeScene: ベンゼンテンプレートの描画_

- assert benzene_button is not None
- assert scene.mode == 'template_benzene'
- assert len(window.data.atoms) == 6
- assert len(window.data.bonds) == 6
- assert orders == [1, 1, 1, 2, 2, 2]

### test_open_settings_dialog
_MainWindow: 設定ダイアログを開くテスト_


### test_toggle_measurement_mode
_MainWindow: 3D測定モードのトグルテスト_

- assert window.measurement_mode == False
- assert measurement_action is not None
- assert window.measurement_mode == True
- assert window.statusBar().currentMessage().startswith('Measurement mode enabled')
- assert window.measurement_mode == False
- assert window.statusBar().currentMessage() == 'Measurement mode disabled.'

### test_toggle_3d_edit_mode
_MainWindow: 3Dドラッグモードのトグルテスト_

- assert window.is_3d_edit_mode == False
- assert edit_3d_action is not None
- assert window.is_3d_edit_mode == True
- assert window.statusBar().currentMessage() == '3D Drag Mode: ON.'
- assert window.is_3d_edit_mode == False
- assert window.statusBar().currentMessage() == '3D Drag Mode: OFF.'

### test_add_remove_hydrogens
_編集: 水素の追加/削除メニュー_

- assert len(window.data.atoms) == 1
- assert len(window.data.atoms) == 5
- assert len(window.data.bonds) == 4
- assert symbols.count('H') == 4
- assert len(window.data.atoms) == 1
- assert len(window.data.bonds) == 0
- assert window.data.atoms[0]['symbol'] == 'C'

### test_2d_cleanup
_2Dクリーンアップ: ボタンクリックで座標が変更されるか_

- assert not (pos0_before == pos0_after and pos1_before == pos1_after)
- assert '2D structure optimization successful' in window.statusBar().currentMessage()

### test_3d_viewer_mode_mol
_3Dビューアモード: 実際のMOLファイル読み込みとUI遷移の統合テスト_

- assert window.current_mol is not None
- assert window.current_mol.GetNumAtoms() == 1
- assert window.is_2d_editable == False
- assert window.cleanup_button.isEnabled() == False
- assert get_button(window.toolbar, 'N (n)').isEnabled() == False
- assert window.optimize_3d_button.isEnabled() == True
- assert window.export_button.isEnabled() == True
- assert window.analysis_action.isEnabled() == True

### test_open_3d_edit_dialogs
_3D編集: 3D編集ダイアログが起動するか_

- assert window.current_mol is not None
- assert window.translation_action.isEnabled() == True
- assert window.align_menu.isEnabled() == True
- assert window.planarize_action.isEnabled() == True
- assert len(window.active_3d_dialogs) == 0

### test_save_project_as
_プロジェクト保存: "Save Project As..." のテスト_

- assert window.has_unsaved_changes == False
- assert window.current_file_path == '/fake/save.pmeprj'
- assert 'Project saved to' in window.statusBar().currentMessage()

### test_open_project
_プロジェクト読み込み: "Open Project..." (.pmeprj) のテスト_

- assert len(window.data.atoms) == 2
- assert len(window.data.bonds) == 1
- assert 0 in window.data.atoms
- assert 1 in window.data.atoms
- assert window.data.atoms[0]['symbol'] == 'C'
- assert (0, 1) in window.data.bonds
- assert window.current_mol is None
- assert 'Project loaded from' in window.statusBar().currentMessage()

### test_toggle_3d_atom_info
_3D原子情報表示: ID, 座標, シンボル表示の切り替えテスト_

- assert window.current_mol is not None
- assert window.atom_info_display_mode == 'id'
- assert window.current_atom_info_labels is not None
- assert window.atom_info_display_mode == 'coords'
- assert window.atom_info_display_mode is None
- assert window.current_atom_info_labels is None

### test_user_template_dialog_save_and_use
_ユーザーテンプレート: ダイアログを開き、現在の構造を保存し、使用するテスト_

- assert action_open_dialog is not None
- assert len(window.data.atoms) == 1
- assert QMessageBox.information.called
- assert any(("Template 'test' saved" in str(a) for a in called_args))
- assert len(window.data.atoms) == 2
- assert 1 in window.data.atoms
- assert window.data.atoms[1]['symbol'] == 'C'
- assert window.data.atoms[1]['pos'].x() > 50

### test_implicit_hydrogens_update
_暗黙の水素: 描画操作後に自動更新されるかのテスト_

- assert len(window.data.atoms) == 1
- assert atom_item.implicit_h_count == 4
- assert len(window.data.atoms) == 2
- assert len(window.data.bonds) == 1
- assert window.data.atoms[0]['item'].implicit_h_count == 3
- assert window.data.atoms[1]['item'].implicit_h_count == 3

### test_drag_drop_mol_file_on_3d_view
_D&D: 3Dビュー領域への .mol ファイルドロップ (モック)_


### test_drag_drop_mol_file_on_2d_view
_D&D: 2Dビュー領域への .mol ファイルドロップ (モック)_


### test_project_save_load_round_trip
_プロジェクト保存/読込: 保存したファイルを実際に読み込んで復元を確認する統合テスト_

- assert len(window.data.atoms) == 2
- assert window.data.atoms[id1]['symbol'] == 'N'
- assert (0, 1) in window.data.bonds
- assert save_file.exists()
- assert window.has_unsaved_changes == False
- assert len(window.data.atoms) == 0
- assert len(window.data.atoms) == 2
- assert len(window.data.bonds) == 1
- assert symbols == ['C', 'N']
- assert len(window.data.bonds) == 1
- assert bond_data['order'] == 1

### test_file_import_smiles_error
_SMILESインポート: 不正なSMILES入力時のエラーハンドリングテスト_

- assert dialog_called or 'Failed' in status_msg or 'Error' in status_msg or ('Invalid' in status_msg)

### test_undo_redo_boundary
_Undo/Redo: スタック境界(空のスタックへの操作)のテスト_

- assert len(window.undo_stack) <= 1
- assert len(window.data.atoms) == 0
- assert len(window.data.atoms) == 0

### test_import_invalid_mol_file
_ファイル読込: 破損したMOLファイルのエラーハンドリング_

- assert mock_msg_box.called or failed_status
- assert len(window.data.atoms) == 0

### test_clear_2d_editor_cancel
_全消去: 確認ダイアログでキャンセルのテスト_

- assert len(window.data.atoms) == 1
- assert window.has_unsaved_changes == True

### test_clipboard_copy_empty_selection
_コピー: 選択なしでのコピー操作の安全性テスト_


## tests/gui/test_plugin_manager.py

### test_init
_Test initialization of PluginManager._

- assert os.path.exists(plugin_manager.plugin_dir)
- assert plugin_manager.plugins == []

### test_install_and_discover_single_file
_Test installing and discovering a single-file plugin._

- assert success
- assert (tmp_path / 'plugins' / 'test_plugin.py').exists()
- assert len(plugins) == 1
- assert p['name'] == 'Test Plugin'
- assert p['version'] == '1.0'
- assert p['description'] == 'A simple test plugin'
- assert p['status'] == 'Loaded'

### test_plugin_registration
_Test that a plugin can register actions via the context._

- assert len(plugin_manager.menu_actions) == 1
- assert action['plugin'] == 'Action Plugin'
- assert action['text'] == 'Test Action'

### test_install_zip
_Test installing a plugin from a ZIP file._

- assert success
- assert extracted_dir.exists()
- assert (extracted_dir / '__init__.py').exists()
- assert len(plugins) == 1
- assert plugins[0]['name'] == 'Zipped Plugin'

### test_ast_metadata_parsing
_Test the safe metadata extraction using AST._

- assert info['name'] == 'AST Plugin'
- assert info['version'] == '1.2.3'
- assert info['author'] == 'Test Author'
- assert info['description'] == 'Docstring description.'
