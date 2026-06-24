# Test Assertions Catalog

## tests/unit/test_about_dialog.py

### test_about_dialog_initialization
_Verify AboutDialog initializes correctly with the injected host._

- assert dialog.windowTitle() == 'About MoleditPy'
- assert dialog.image_label is not None

### test_about_dialog_easter_egg
_Verify the easter egg triggers correctly on right click._

- mock_parser_host.edit_actions_manager.clear_all.assert_called_once()
- mock_parser_host.string_importer_manager.load_from_smiles.assert_called_once_with('C1=CN=C(N=C1)C2=NC=CC=N2')

### test_about_dialog_ignore_left_click
_Verify left clicks do not trigger the easter egg._

- mock_parser_host.clear_all.assert_not_called()
- event.ignore.assert_called_once()

## tests/unit/test_align_plane_dialog.py

### TestApplyPlaneAlignGuard.test_fewer_than_three_atoms_shows_warning
_Warning shown when fewer than 3 atoms are selected for plane alignment._

- mb.warning.assert_called_once()

### TestApplyPlaneAlignGuard.test_zero_atoms_shows_warning
_Warning shown when no atoms are selected for plane alignment._

- mb.warning.assert_called_once()

### TestApplyPlaneAlignMath.test_xy_align_reduces_z_variance
_XY alignment brings all selected atom Z-coordinates close to zero._

- assert abs(after[i][2]) < 0.5

### TestApplyPlaneAlignMath.test_apply_calls_draw_molecule_3d
_Plane alignment triggers a 3D redraw after modifying geometry._

- mw.view_3d_manager.draw_molecule_3d.assert_called()

### TestApplyPlaneAlignMath.test_align_plane_with_move_to_zero_plane_true
_move_to_zero_plane=True translates selected atoms to Z=0 after alignment._

- assert pos[i][2] == pytest.approx(0.0, abs=1e-05)

### TestApplyPlaneAlignMath.test_align_plane_with_move_to_zero_plane_false
_move_to_zero_plane=False keeps the centroid Z offset after alignment._

- assert abs(pos[i][2]) > 1.0

### TestApplyPlaneAlignMath.test_align_plane_already_aligned_with_move_to_zero_plane_true
_Already-aligned plane with move_to_zero still zeroes the Z offset._

- assert pos[i][2] == pytest.approx(0.0, abs=1e-05)

## tests/unit/test_alignment_dialog.py

### TestOnAtomPicked.test_first_pick_adds_atom
_Picking an atom for the first time adds it to selected_atoms._

- assert 0 in dlg.selected_atoms

### TestOnAtomPicked.test_second_pick_adds_atom
_Picking a second distinct atom accumulates both atoms in selected_atoms._

- assert dlg.selected_atoms == {0, 1}

### TestOnAtomPicked.test_third_pick_capped_at_two
_A third pick is ignored; AlignmentDialog accepts at most two atoms._

- assert len(dlg.selected_atoms) == 2
- assert 2 not in dlg.selected_atoms

### TestOnAtomPicked.test_repick_deselects
_Picking the same atom twice removes it from selected_atoms._

- assert 0 not in dlg.selected_atoms

### TestOnAtomPicked.test_enables_apply_when_two_atoms
_Selecting two atoms enables the apply button._

- assert dlg.apply_button.isEnabled()

### TestPreselectedAtoms.test_preselected_loaded
_Preselected atoms are loaded into selected_atoms on construction._

- assert dlg.selected_atoms == {0, 1}

### TestPreselectedAtoms.test_preselected_capped_at_two
_More than two preselected atoms are silently capped to two._

- assert len(dlg.selected_atoms) == 2

### TestClearSelection.test_clears_atoms_and_disables_apply
_clear_selection empties selected_atoms and disables the apply button._

- assert len(dlg.selected_atoms) == 0
- assert not dlg.apply_button.isEnabled()

### TestUpdateDisplay.test_zero_atoms_label
_update_display with no selection disables the apply button._

- assert not dlg.apply_button.isEnabled()

### TestUpdateDisplay.test_one_atom_label_contains_symbol
_update_display with one atom shows the element symbol and disables apply._

- assert sym in dlg.selection_label.text()
- assert not dlg.apply_button.isEnabled()

### TestUpdateDisplay.test_two_atoms_enables_apply
_update_display with two atoms enables apply and shows count in label._

- assert dlg.apply_button.isEnabled()
- assert '2' in dlg.selection_label.text()

### TestApplyAlignmentGuard.test_fewer_than_two_atoms_shows_warning
_apply_alignment shows a warning when fewer than two atoms are selected._

- mb.warning.assert_called_once()

### TestApplyAlignmentGuard.test_zero_atoms_shows_warning
_apply_alignment with empty selection shows a warning._

- mb.warning.assert_called_once()

### TestApplyAlignmentMath.test_x_axis_alignment_atom1_at_origin
_After X-alignment, atom1 must be at origin._

- assert pos[0] == pytest.approx([0.0, 0.0, 0.0], abs=1e-05)

### TestApplyAlignmentMath.test_x_axis_alignment_atom2_on_x_axis
_After X-alignment, the molecule is aligned parallel to X-axis, rotated about centroid._

- assert pos[1][1] == pytest.approx(1.5, abs=1e-05)
- assert pos[1][2] == pytest.approx(0.0, abs=1e-05)
- assert pos[1][0] > 0

### TestApplyAlignmentMath.test_z_axis_alignment_atom2_on_z_axis
_After Z-alignment, the molecule is aligned parallel to Z-axis, rotated about centroid._

- assert pos[1][0] == pytest.approx(1.5, abs=1e-05)
- assert pos[1][1] == pytest.approx(0.0, abs=1e-05)
- assert pos[1][2] > 0

### TestApplyAlignmentMath.test_apply_pushes_undo
_apply_alignment calls push_undo_state to record the geometry change._

- mw.edit_actions_manager.push_undo_state.assert_called()

### TestApplyAlignmentMath.test_alignment_with_move_to_origin_true
_When move_to_origin is True, the first atom is moved to the origin._

- assert pos[0] == pytest.approx([0.0, 0.0, 0.0], abs=1e-05)
- assert pos[1] == pytest.approx([3.0, 0.0, 0.0], abs=1e-05)

### TestApplyAlignmentMath.test_alignment_with_move_to_origin_false
_When move_to_origin is False, the molecule rotates about its centroid (centroid preserved)._

- assert pos[0] == pytest.approx([-0.5, 3.5, 3.0], abs=1e-05)
- assert pos[1] == pytest.approx([2.5, 3.5, 3.0], abs=1e-05)

### TestApplyAlignmentMath.test_already_aligned_with_move_to_origin_true
_If already aligned, setting move_to_origin to True should still translate the first atom to origin._

- assert pos[0] == pytest.approx([0.0, 0.0, 0.0], abs=1e-05)

## tests/unit/test_alt_template_bypass.py

### test_update_template_preview_allows_snapping_but_bypasses_fusing_when_alt_pressed
_Alt key allows cursor snapping but suppresses template fusing in preview._

- assert len(scene.find_atom_near_args) == 1
- assert passed_points[1].x() == 10.0
- assert passed_points[1].y() == 20.0

### test_add_molecule_fragment_bypasses_fusing_when_alt_pressed
_Alt key bypasses fusing so all template fragment atoms are created fresh._

- assert scene.create_atom.call_count == 6

## tests/unit/test_app_logic.py

### test_ez_preservation_logic
_Verify worker preserves explicit E/Z labels even if RDKit might lose them._

- assert bond.GetStereo() == Chem.BondStereo.STEREOZ

### test_molecular_data_fallback_serialization
_Verify MolecularData uses manual string construction if RDKit fails (fallback logic)._

- assert 'V2000' in counts_line
- assert 'MoleditPy' in mol_block
- assert 'X  ' in mol_block
- assert 'C  ' in mol_block

### test_coordinate_mapping_primary
_Verify primary RDKit conversion logic uses 'pos' attribute correctly._

- assert len(atom_lines) == 2
- assert abs(x_coord - expected_x) < 0.0001

### test_coordinate_mapping_fallback
_Verify fallback serialization (reverted logic) uses atom['item'].pos()._

- assert 'MoleditPy' in mol_block
- assert len(atom_lines) == 2
- assert abs(x_coord - expected_x) < 0.0001

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

- assert isinstance(mol_3d_binary, bytes) and len(mol_3d_binary) > 0
- assert restored.GetNumAtoms() == mol.GetNumAtoms()
- assert 'mol_3d_atom_ids' in state

### test_get_current_state_includes_version
_State should always include version string._

- assert 'version' in state

### test_get_current_state_captures_constraints
_3D constraints should be serialized in JSON-safe format._

- assert len(state['constraints_3d']) == 2
- assert isinstance(c, list)

### test_get_current_state_empty
_Empty state should produce valid structure with no atoms/bonds._

- assert state['atoms'] == {}
- assert state['bonds'] == {}

### test_json_roundtrip_preserves_atoms
_Atoms should survive JSON serialization round-trip._

- assert len(mock_parser_host.data.atoms) == 2
- assert symbols == {'C', 'N'}
- assert len(charged) == 1
- assert charged[0]['charge'] == 1

### test_json_roundtrip_preserves_bonds
_Bond order should survive JSON round-trip._

- assert len(mock_parser_host.data.bonds) == 1
- assert bond['order'] == 2

### test_json_roundtrip_preserves_radical
_Radical electrons should survive JSON round-trip._

- assert len(mock_parser_host.data.atoms) == 1
- assert next(iter(mock_parser_host.data.atoms.values()))['radical'] == 2

### test_pmeprj_serialization_roundtrip
_Test full project serialization/deserialization (PMEPRJ)._

- assert json_data['format'] == 'PME Project'
- assert '2d_structure' in json_data
- assert '3d_structure' in json_data
- assert json_data['3d_structure']['num_conformers'] == 1
- assert len(mw.state_manager.data.atoms) == 2
- assert mw.state_manager.data.atoms[aid1]['symbol'] == 'C'
- assert mw.state_manager.data.atoms[aid2]['symbol'] == 'O'
- assert mw.state_manager.data.atoms[aid1]['pos'][0] == 10
- assert mw.view_3d_manager.current_mol is not None
- assert mw.view_3d_manager.current_mol.GetNumAtoms() == mol.GetNumAtoms()
- assert mw.view_3d_manager.current_mol.GetAtomWithIdx(0).GetIntProp('_original_atom_id') == aid1
- assert mw.view_3d_manager.current_mol.GetAtomWithIdx(1).GetIntProp('_original_atom_id') == aid2
- assert mw.view_3d_manager.atom_positions_3d is not None
- assert mw.view_3d_manager.atom_positions_3d.shape == (mol.GetNumAtoms(), 3)
- assert len(mw.edit_3d_manager.constraints_3d) == 1
- assert mw.edit_3d_manager.constraints_3d[0][0] == 'DISTANCE'
- assert mw.edit_3d_manager.constraints_3d[0][2] == 1.43
- assert mw.compute_manager.last_successful_optimization_method == 'MMFF94s'
- assert mw._preserved_plugin_data['TestPlugin']['val'] == 42

### test_undo_state_binary_roundtrip
_Test the internal binary state serialization used for Undo/Redo._

- assert 'mol_3d' in state
- assert isinstance(state['mol_3d'], bytes)
- assert state['mol_3d_atom_ids'][0] == aid
- assert mw.state_manager.data.atoms[aid]['symbol'] == 'N'
- assert mw.view_3d_manager.current_mol.GetAtomWithIdx(0).GetIntProp('_original_atom_id') == aid

### test_legacy_version_handling
_Verify that version mismatch warnings are triggered (but don't crash)._


## tests/unit/test_atom_bond_items.py

### TestAtomItem.test_init
_Verify AtomItem initialization with ID, symbol and position._

- assert atom_item.symbol == 'C'
- assert atom_item.atom_id == 1
- assert atom_item.pos().x() == 0.0
- assert atom_item.pos().y() == 0.0
- assert atom_item.flags() & atom_item.GraphicsItemFlag.ItemIsMovable
- assert atom_item.flags() & atom_item.GraphicsItemFlag.ItemIsSelectable

### TestAtomItem.test_paint_mock
_Test paint logic by mocking QPainter_

- assert mock_painter.drawText.called
- assert text_drawn

### TestAtomItem.test_paint_radical
_Test painting logic for radicals_

- assert mock_painter.drawEllipse.called
- assert mock_painter.drawEllipse.call_count >= 2

### TestAtomItem.test_paint_selection_hover
_Test painting logic for selection and hover highlights_

- assert mock_painter.drawRect.called
- assert mock_painter.drawRect.called
- assert mock_painter.drawRect.called

### TestAtomItem.test_paint_implicit_hydrogens_and_flipping
_Test painting logic for implicit hydrogens and text alignment_

- assert found_subscript
- assert found_flipped_text

### TestBondItem.test_init
_Verify BondItem initialization with atom partners and default order._

- assert bond_item.atom1.atom_id == 1
- assert bond_item.atom2.atom_id == 2
- assert bond_item.order == 1

### TestBondItem.test_update_position
_Verify that the bond line updates correctly when atom positions change._

- assert line.p1() == QPointF(0.0, 0.0)
- assert line.p2() == QPointF(5.0, 5.0)

### TestBondItem.test_geometric_shortening
_Verify that bonds are geometrically shortened to not cross atom labels._

- assert line.p1().x() == pytest.approx(5.0, abs=0.1)
- assert line.p2().x() == pytest.approx(15.0, abs=0.1)
- assert line.p1().y() == pytest.approx(0.0, abs=0.1)
- assert line.p2().y() == pytest.approx(0.0, abs=0.1)

### TestBondItem.test_set_bond_order
_Verify that the bond order can be changed and is reflected in the item state._

- assert bond_item.order == 2
- assert bond_item.order == 3

### TestBondItem.test_paint_mock
_Test paint logic for bond_

- assert mock_painter.drawLine.called
- assert mock_painter.drawLine.call_count >= 1

### TestBondItem.test_paint_ring_double_bond
_Test painting logic for double bond in a ring_

- assert mock_painter.drawLine.call_count >= 2

### TestBondItem.test_ring_inner_double_bond_shortening_preserves_benzene
_Benzene-like 60 degree ring bonds keep the established inner-line length._

- assert factor == pytest.approx(0.8)

### TestBondItem.test_ring_inner_double_bond_shortening_handles_small_rings
_Small rings with sharper spans shorten the inner double bond more._

- assert factor == pytest.approx(0.55)

### TestBondItem.test_bounding_rect_ez_label
_Test boundingRect expansion for E/Z labels_

- assert rect_stereo.width() >= rect_no_stereo.width()
- assert rect_stereo.height() >= rect_no_stereo.height()

## tests/unit/test_atom_picking.py

### test_pick_atom_index_from_screen_hits_projected_atom_edge
_Click on a projected atom edge returns that atom's index._

- assert pick_atom_index_from_screen(_view(), (111, 100), _Mol()) == 0

### test_pick_atom_index_from_screen_returns_none_for_background
_Click on empty background returns None._

- assert pick_atom_index_from_screen(_view(), (200, 200), _Mol()) is None

### test_pick_atom_index_from_screen_vectorized_success
_Vectorized implementation returns the correct atom index on a hit._

- assert pick_atom_index_from_screen_vectorized(_view(), (111, 100), _Mol()) == 0

### test_pick_atom_index_from_screen_sequential_success
_Sequential implementation returns the correct atom index on a hit._

- assert pick_atom_index_from_screen_sequential(_view(), (111, 100), _Mol()) == 0

### test_pick_atom_index_from_screen_fallback
_Falls back to sequential picking when the camera is unavailable._

- assert pick_atom_index_from_screen(view_obj, (111, 100), _Mol()) == 0

## tests/unit/test_atom_placement_logic.py

### test_placement_0_neighbors
_Zero bonds: should place default (up)._

- assert offset.x() == pytest.approx(0)
- assert offset.y() == pytest.approx(-L)

### test_placement_1_neighbor
_One bond: should rotate 60 degrees clockwise from existing bond vector._

- assert offset.x() == pytest.approx(10.0)
- assert offset.y() == pytest.approx(17.32, abs=0.01)

### test_placement_3_neighbors_balanced
_Three balanced neighbors: sum is near zero, should use fallback (45 deg offset)._

- assert offset.x() == pytest.approx(L * 0.7071)
- assert offset.y() == pytest.approx(-L * 0.7071)

### test_placement_3_neighbors_unbalanced
_Three unbalanced neighbors: should place in direction opposite to the sum._

- assert offset.x() == pytest.approx(0)
- assert offset.y() == pytest.approx(-L)

### test_placement_2_neighbors_skeleton
_Two neighbors: should continue skeleton (opposite to average bond vector)._

- assert offset.x() == pytest.approx(0)
- assert offset.y() == pytest.approx(-L)

### test_placement_1_neighbor_anticlockwise
_One bond with placement_direction_clockwise=False: rotate -60 degrees (anticlockwise)._

- assert offset.x() == pytest.approx(10.0)
- assert offset.y() == pytest.approx(-17.32, abs=0.01)

### test_placement_alkyne_straight_continuation_target_order_3
_If target_order is 3 (triple bond), the angle should be 0 (straight continuation)._

- assert offset.x() == pytest.approx(L)
- assert offset.y() == pytest.approx(0.0, abs=0.01)

### test_placement_alkyne_straight_continuation_existing_bond_order_3
_If the existing bond order is 3, the angle should be 0 (straight continuation)._

- assert offset.x() == pytest.approx(L)
- assert offset.y() == pytest.approx(0.0, abs=0.01)

## tests/unit/test_base_picking_dialog.py

### test_key_enter_clicks_apply_button_if_enabled
_Pressing Enter clicks the apply button when it is enabled._

- apply_btn.click.assert_called_once()
- ev.accept.assert_called_once()

### test_key_enter_does_not_click_disabled_apply_button
_Pressing Enter does not click the apply button when it is disabled._

- apply_btn.click.assert_not_called()

### test_key_enter_no_apply_button_does_not_raise
_Pressing Enter when apply_button is absent does not raise._


### test_key_none_event_returns_early
_None event passed to keyPressEvent returns early without raising._


### test_key_other_key_passes_to_super
_Non-Enter key events are forwarded to QDialog.keyPressEvent without error._


### test_cleanup_clears_labels_and_disables_picking
_closeEvent, reject, and accept all clear atom labels and disable picking._

- mock_clear.assert_called_once()
- mock_disable.assert_called_once()

### test_update_molecule_geometry_array_updates_conformer
__update_molecule_geometry with an ndarray updates the conformer positions._

- assert updated[0] == pytest.approx([10.0, 20.0, 30.0], abs=0.0001)
- assert dlg._molecule_modified is True

### test_update_molecule_geometry_dict_form
__update_molecule_geometry accepts a dict mapping atom index to position._

- assert updated[0] == pytest.approx([5.0, 6.0, 7.0], abs=0.0001)

### test_update_molecule_geometry_calls_draw_molecule_3d
__update_molecule_geometry triggers a 3D redraw after updating positions._

- mw.view_3d_manager.draw_molecule_3d.assert_called_once_with(mol)

### test_update_molecule_geometry_updates_cache
__update_molecule_geometry also updates the atom_positions_3d cache._

- assert mw.view_3d_manager.atom_positions_3d[1] == pytest.approx([99.0, 0.0, 0.0], abs=0.0001)

### test_push_undo_calls_push_undo_state
__push_undo calls push_undo_state and clears the _molecule_modified flag._

- mw.edit_actions_manager.push_undo_state.assert_called_once()
- assert dlg._molecule_modified is False

### test_push_undo_no_state_manager_no_crash
__push_undo does not raise when state_manager is absent._


### test_done_pushes_undo_if_molecule_modified
_done() calls _push_undo when _molecule_modified is True._

- mock_undo.assert_called_once()

### test_done_skips_undo_if_not_modified
_done() skips _push_undo when _molecule_modified is False._

- mock_undo.assert_not_called()

## tests/unit/test_benzene_placement_shortcut.py

### test_benzene_shortcut_on_atom
_Test one-shot benzene placement when cursor is over an atom._

- scene._calculate_polygon_from_edge.assert_called_once()
- scene.add_molecule_fragment.assert_called_once()
- scene.window.edit_actions_manager.push_undo_state.assert_called_once()
- event.accept.assert_called_once()

### test_benzene_shortcut_on_bond
_Test one-shot benzene placement when cursor is over a bond._

- assert called_kwargs.get('use_existing_length') is True
- scene.add_molecule_fragment.assert_called_once()

### test_benzene_shortcut_empty_space
_Test mode switch when cursor is over empty space._

- scene.window.ui_manager.set_mode_and_update_toolbar.assert_called_with('template_benzene')
- event.accept.assert_called_once()
- scene.add_molecule_fragment.assert_not_called()

## tests/unit/test_benzene_rotation.py

### test_calculate_6ring_rotation_empty
_Test with no existing bonds._

- assert rot == 0

### test_calculate_6ring_rotation_single_edge_single
_Test fusing on a single bond (order 1). Should prefer alternating (template double)._

- assert rot % 2 == 0

### test_calculate_6ring_rotation_single_edge_double
_Test fusing on a double bond. Should prefer 1 or 2 (alternating or matching)._

- assert rot % 2 == 0

### test_calculate_6ring_rotation_multi_edge_fused
_Test fusing on two adjacent edges (naphthalene-like)._

- assert rot == 0

### test_calculate_6ring_rotation_connection_safety
_Verify that 'safe connection' scoring prioritizes rotations where template single bonds connect to fusion points._

- assert rot % 2 == 0

## tests/unit/test_calculation_worker_core.py

### test_calculation_worker_init
_CalculationWorker initialises as a QObject with halt attributes unset._

- assert isinstance(worker, QObject)
- assert getattr(worker, 'halt_ids', None) is None
- assert not getattr(worker, 'halt_all', False)

### test_calculation_worker_explicit_stereo_m_cfg
_M CFG stereo block in V2000 molfile is preserved through direct conversion._

- assert len(finish_captor.emitted_values) > 0
- assert bond.GetStereo() == Chem.BondStereo.STEREOE

### test_calculation_worker_error_empty_input
_Empty mol block emits an error signal indicating no atoms to convert._

- assert any(('No atoms to convert' in str(val) for val in error_captor.emitted_values))

### test_calculation_worker_none_conversion_mode_defaults_to_fallback
_conversion_mode=None falls through to the RDKit workflow path._

- assert mock_rdkit.called

### test_calculation_worker_safe_helpers_halted
_halt_all=True prevents the finished signal from being emitted._

- assert len(finish_captor.emitted_values) == 0

### test_calculation_worker_rdkit_embedding_fail_fallback
_Embedding failure is reported via status/error signal and does not raise._

- assert any((keyword in all_msgs_str for keyword in ['embedding failed', 'conversion failed', 'bounds fail']))

### test_optimize_only_mmff94s
_optimize_only with MMFF94s sets the _pme_optimization_method property._

- assert len(finish_captor.emitted_values) > 0
- assert res_mol.HasProp('_pme_optimization_method')
- assert res_mol.GetProp('_pme_optimization_method').upper() == 'MMFF94S_RDKIT'

### test_optimize_only_uff
_optimize_only with UFF_RDKIT sets the _pme_optimization_method property._

- assert len(finish_captor.emitted_values) > 0
- assert res_mol.GetProp('_pme_optimization_method') == 'UFF_RDKIT'

### test_collision_avoidance_trigger
_Overlapping atoms in direct mode trigger collision avoidance separation._

- assert len(finish_captor.emitted_values) > 0
- assert dist > 0.01

### test_iterative_optimize_halt
__iterative_optimize raises WorkerHaltError when halt_ids contains the worker id._


### test_obabel_optimization_flow
_optimize_only with UFF_OBABEL calls the obabel optimization helper._

- assert mock_opt.called
- assert len(finish_captor.emitted_values) > 0
- assert res_mol.GetProp('_pme_optimization_method') == 'UFF_OBABEL'

### test_intermolecular_interaction_toggle
_Enabling intermolecular interaction optimization draws separate molecules together._

- assert abs(dist_off - 6.0) < 1e-05
- assert dist_on < 5.0

### test_intermolecular_interaction_uff
_UFF with intermolecular interaction disabled keeps molecules at their initial separation._

- assert abs(dist_off - 6.0) < 1e-05

### test_iterative_optimize_mmff_unsupported_returns_false
__iterative_optimize returns False (no UFF fallback) when MMFF props are None._

- assert result is False
- assert any(('Failed to setup' in m for m in status_msgs))

### test_optimize_only_mmff_unsupported_emits_error
_optimize_only with unsupported MMFF (props=None) emits an error signal, not success._

- assert len(finish_captor.emitted_values) == 0
- assert len(error_captor.emitted_values) > 0
- assert 'MMFF' in err_msg.upper() or 'failed' in err_msg.lower()

## tests/unit/test_color_settings_dialog.py

### test_color_settings_dialog_initialization
_Verify ColorSettingsDialog initializes and parses settings correctly._

- assert dialog.windowTitle() == 'CPK Colors'
- assert 'C' in dialog.element_buttons
- assert dialog.changed_cpk == {}

### test_color_settings_dialog_pick_color
_Verify that clicking an element button updates the changed_cpk dictionary._

- assert dialog.changed_cpk['C'] == '#ff0000'
- assert 'background-color: #ff0000' in btn.styleSheet()

### test_color_settings_dialog_reset_all
_Verify reset_all clears changes and sets the reset flag._

- assert dialog._reset_all_flag is True
- assert dialog.changed_cpk == {}

### test_color_settings_dialog_apply_changes
_Verify that apply_changes pushes updates to the parent window settings._

- assert mock_parser_host.settings_dirty is True
- assert mock_parser_host.settings['cpk_colors']['O'] == '#00ff00'
- assert mock_parser_host.settings['ball_stick_bond_color'] == '#112233'

### test_pick_bs_bond_color_valid_updates_changed_bs
_A valid colour from QColorDialog is stored in changed_bs_color._

- assert dialog.changed_bs_color == '#aabbcc'
- assert 'background-color: #aabbcc' in dialog.bs_button.styleSheet()

### test_pick_bs_bond_color_invalid_no_change
_An invalid colour from QColorDialog leaves changed_bs_color as None._

- assert dialog.changed_bs_color is None

### test_reset_all_with_parent_uses_default_bond_color
_reset_all uses the parent's default_settings bond colour when available._

- assert dialog.changed_bs_color == '#123456'
- assert 'background-color: #123456' in dialog.bs_button.styleSheet()

### test_reset_all_without_parent_defaults_to_gray
_reset_all falls back to #7F7F7F when no parent window is set._

- assert dialog.changed_bs_color == '#7F7F7F'
- assert 'background-color: #7f7f7f' in dialog.bs_button.styleSheet().lower()

### test_reset_all_restores_element_buttons
_reset_all clears changed_cpk and sets the reset flag._

- assert dialog.changed_cpk == {}
- assert dialog._reset_all_flag is True

### test_apply_changes_no_parent_returns_early
_apply_changes is a no-op when parent_window is None._


### test_apply_changes_reset_flag_deletes_cpk_colors
_apply_changes with _reset_all_flag removes cpk_colors from settings._

- assert 'cpk_colors' not in parent.init_manager.settings

### test_apply_changes_with_mol_calls_draw_molecule_3d
_apply_changes triggers a 3D redraw when a molecule is loaded._

- parent.view_3d_manager.draw_molecule_3d.assert_called()

### test_apply_changes_reset_flag_resets_bond_color
_apply_changes with _reset_all_flag resets ball-stick bond colour to default._

- assert parent.init_manager.settings.get('ball_stick_bond_color') == '#7F7F7F'

### test_apply_changes_updates_2d_scene_items
_apply_changes calls update_style on 2D scene items that support it._

- item_with_style.update_style.assert_called()

### test_apply_changes_calls_update_cpk_colors
_apply_changes always calls update_cpk_colors_from_settings on the parent._

- parent.init_manager.update_cpk_colors_from_settings.assert_called()

### test_apply_changes_calls_save_settings_when_cpk_changed
_apply_changes calls save_settings when CPK colours have changed._

- parent.init_manager.save_settings.assert_called()

### test_accept_calls_apply_then_super
_accept() calls apply_changes before delegating to QDialog.accept._

- mock_apply.assert_called_once()
- mock_super.assert_called_once()

### test_on_element_clicked_invalid_color_no_change
_on_element_clicked does not update changed_cpk when the colour dialog is cancelled._

- assert 'C' not in dialog.changed_cpk

### test_init_cpk_override_applied_to_button
_CPK colour overrides in current_settings are reflected in button styles on init._

- assert '#112233' in style

## tests/unit/test_compute_logic.py

### test_on_calculation_error_stale
_Test on_calculation_error when the worker is stale (not in active set)._

- compute.statusBar().showMessage.assert_called()
- assert 'stale' in msg.lower() or 'Ignored' in msg

### test_on_calculation_error_basic
_Test on_calculation_error for an ACTIVE worker._

- assert compute.statusBar().showMessage.called
- assert 'Real Error' in compute.statusBar().showMessage.call_args[0][0]

### test_compute_set_optimization_method
_Verify that setting the optimization method updates both settings and internal state._

- assert compute.host.init_manager.settings['optimization_method'] == 'GAFF_OBABEL'
- assert compute.statusBar().showMessage.called
- assert 'Optimization' in msg or 'GAFF_OBABEL' in msg
- assert compute.host.init_manager.optimization_method == 'GAFF_OBABEL'

### test_compute_halt_logic
_Verify that halt_conversion correctly marks active workers for termination._

- assert 'test_id' in compute.halt_ids
- assert len(compute.active_worker_ids) == 0
- assert compute.statusBar().showMessage.called

### test_on_calculation_finished_basic
_Verify that on_calculation_finished correctly processes a finished worker result._

- assert compute.host.view_3d_manager.current_mol == mol
- assert worker_id not in compute.active_worker_ids

### test_check_chemistry_problems_fallback_detects
_Verify that the manual valence fallback correctly identifies overvalent atoms._

- assert c_item.has_problem is True
- assert compute.statusBar().showMessage.called

### test_trigger_conversion_empty
_Verify that trigger_conversion handles empty molecular data gracefully._

- assert compute.statusBar().showMessage.called

### test_trigger_conversion_with_atoms
_Verify that trigger_conversion correctly starts the calculation thread for a valid molecule._

- assert compute.statusBar().showMessage.called

### test_optimize_3d_structure_logic
_Verify the high-level logic of triggering 3D optimization on the current molecule._

- assert compute.statusBar().showMessage.called

### test_on_calculation_finished_worker_id_mismatch
_Verify that on_calculation_finished ignores results from stale or mismatched workers._

- assert compute.host.view_3d_manager.current_mol is None

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

### test_optimize_3d_temp_method_override
_Test optimize_3d_structure with temporary optimization method override._

- assert MockWorker.called
- compute.statusBar().showMessage.assert_any_call('Optimizing 3D structure (MMFF_RDKIT)...')

### test_optimize_3d_mmff94s_success
_Test MMFF94s optimization succeeds._

- assert any(('Optimizing' in msg for msg in compute.get_status_messages()))

### test_optimize_3d_uff_success
_Test UFF optimization succeeds._

- assert any(('Optimizing' in msg for msg in compute.get_status_messages()))

### test_optimize_3d_no_conformer
_Test optimize_3d_structure when molecule has no conformer._

- assert any(('No conformer found' in msg for msg in messages))

### test_optimize_3d_mmff_exception_handling
_Test grace during MMFF exception._

- assert MockWorker.called

### test_optimize_3d_uff_exception_handling
_Test grace during UFF exception._

- assert MockWorker.called

### test_optimize_3d_plugin_method
_Test plugin optimization method — worker should be started (not rejected)._

- assert not any(('not available' in msg for msg in msgs))
- assert any(('Optimizing' in msg for msg in msgs))

### test_optimize_3d_plugin_failure
_Test unavailable/unknown method shows correct error._

- assert any(('not available' in msg for msg in msgs))

### test_optimize_3d_mmff_fallback_success
_Test MMFF fallback to ForceField API when basic optimization fails._

- assert any(('Optimizing' in msg for msg in compute.get_status_messages()))

### test_optimize_3d_uff_fallback_failure
_Test UFF fallback failure handles gracefully._

- assert MockWorker.called

### test_optimize_3d_unavailable_method
_Test error when optimization method is unavailable._

- assert any(('Selected optimization' in msg for msg in messages))

### test_on_calculation_error_mmff_fail_shows_uff_dialog
_When MMFF fails, on_calculation_error shows a UFF retry dialog instead of auto-fallback._

- mock_dialog.assert_called_once()
- assert 'UFF' in call_args[2]

### test_on_calculation_error_mmff_fail_user_accepts_uff
_When user accepts UFF retry after MMFF failure, optimize_3d_structure is called with UFF_RDKIT._

- mock_optimize.assert_called_once_with('UFF_RDKIT')

### test_on_calculation_error_mmff_fail_user_declines_uff
_When user declines UFF retry, a critical error dialog is shown instead._

- mock_critical.assert_called_once()

### test_on_calculation_error_non_mmff_no_uff_dialog
_Non-MMFF errors do not trigger the UFF retry dialog._

- mock_dialog.assert_not_called()

### test_on_calculation_finished_collision_single_frag
_Test collision logic is skipped for single fragment._

- assert not mock_adjust.called
- assert not any(('Detecting collisions' in msg for msg in compute.get_status_messages()))

### test_molecular_data_radical_transfer
_Test that radical electrons are transferred correctly to RDKit mol._

- assert mol is not None
- assert mol.GetAtomWithIdx(0).GetNumRadicalElectrons() == 1

### test_app_state_radical_and_constraint_preservation
_Test that radicals and constraints are preserved through state round-trip._

- assert state['atoms'][0]['radical'] == 2
- assert state['constraints_3d'][0] == ['DISTANCE', [0, 1], 1.5, 100000.0]
- assert 0 in compute.data.atoms
- assert compute.data.atoms[0]['radical'] == 2
- assert compute.constraints_3d == [('DISTANCE', (0, 1), 1.5, 100000.0)]

### test_app_state_original_atom_id_preservation
_Test that _original_atom_id is preserved in 3D molecule state round-trip._

- assert 123 in state['mol_3d_atom_ids']
- assert compute.current_mol is not None
- assert compute.current_mol.GetAtomWithIdx(0).HasProp('_original_atom_id')
- assert compute.current_mol.GetAtomWithIdx(0).GetIntProp('_original_atom_id') == 123

### test_optimize_3d_method_persistence
_Test that the optimization method is recorded after success._

- assert compute.last_successful_optimization_method == 'MMFF94s (RDKit)'

### test_trigger_conversion_early_exits
_Test early exits in trigger_conversion (empty mol, etc.)._

- assert compute.host.view_3d_manager.current_mol is None
- assert any(('3D view cleared' in msg for msg in compute.get_status_messages()))
- assert mock_fallback.called

### test_check_chemistry_problems_fallback
_Test the manual valence check when RDKit fails._

- assert mock_item.has_problem == True
- assert any(('chemistry problems found' in msg for msg in msgs))

### test_trigger_conversion_happy_path
_Test trigger_conversion follows the success path until worker setup._

- assert any(('Calculating 3D structure' in msg for msg in msgs))

### test_trigger_conversion_stereo_enhancement
_Test trigger_conversion stereo enhancement logic for E/Z bonds._

- assert mock_timer.called

### test_halt_conversion
_Test halt_conversion clears active workers and restores UI._

- assert compute.active_worker_ids == set()
- assert 1 in compute.halt_ids
- assert compute.host.init_manager.convert_button.setText.called

### test_on_calculation_error_updated
_Test on_calculation_error with correct message formatting._

- assert any(('Error: Test Error' in msg for msg in msgs))
- assert worker_id not in compute.active_worker_ids

### test_trigger_conversion_chemistry_problem_detection
_Test trigger_conversion detects and flags chemistry problems (valence)._

- assert any(('chemistry problems found' in msg for msg in msgs))
- assert mock_item.has_problem == True

### test_trigger_conversion_fragment_message_exact
_Test verification of the exact status bar message for multiple fragments._

- assert any(('Converting 3 molecules' in msg for msg in all_messages))

### test_trigger_conversion_to_mol_block_priority
_Test that data.to_mol_block() is used preferentially over RDKit's generation._

- assert mock_to_block.called
- assert not mock_rdkit_block.called
- assert len(args) >= 2
- assert call_args is not None
- assert sent_block == custom_block

### test_trigger_conversion_uses_default_when_temp_mode_is_none
_Preinitialized temp conversion mode should not suppress the saved setting._

- assert options['conversion_mode'] == 'fallback'

### test_trigger_conversion_ez_stereo_injection
_Test that M CFG lines are injected for E/Z stereo bonds._

- assert 'M  CFG' in sent_block
- assert 'M  CFG  1   2   2' in sent_block

### test_on_calculation_error_uff_fallback_temporary
_Verify that UFF fallback uses _temp_optimization_method and doesn't change persistent setting._

- mock_optimize.assert_called_once_with('UFF_RDKIT')
- assert compute.optimization_method == 'MMFF_RDKIT'

## tests/unit/test_constants.py

### test_constants_version_from_metadata
_Test when importlib.metadata successfully finds the package version._

- assert constants.VERSION == '1.2.3'
- assert mock_version.call_count >= 1

### test_constants_version_metadata_package_not_found
_Test fallback to pyproject.toml when PackageNotFoundError is raised._

- assert constants.VERSION == '2.3.4'

### test_constants_version_import_error_fallback
_Test fallback to pyproject.toml when importlib.metadata cannot be imported (ImportError)._

- assert constants.VERSION == '3.4.5'

### test_constants_version_all_fail_returns_unknown
_Test that it falls back to 'Unknown' when all lookup methods fail._

- assert constants.VERSION == 'Unknown'

### test_constants_version_file_exception_returns_unknown
_Test that any reading file exception is handled and returns 'Unknown'._

- assert constants.VERSION == 'Unknown'

### test_restore_state
_Helper to restore constants state after reloading it under mocks._


## tests/unit/test_constrained_optimization_dialog.py

### TestAtomSelection.test_pick_adds_atom
_Clicking an atom appends its index to selected_atoms._

- assert dlg.selected_atoms == [3]

### TestAtomSelection.test_pick_same_atom_toggles_off
_Clicking the same atom twice removes it from the selection._

- assert dlg.selected_atoms == []

### TestAtomSelection.test_pick_four_atoms
_Picking four distinct atoms fills selected_atoms to capacity._

- assert dlg.selected_atoms == [0, 1, 2, 3]

### TestAtomSelection.test_fifth_atom_evicts_first_fifo
_Adding a 5th atom should drop the oldest (index 0)._

- assert dlg.selected_atoms == [1, 2, 3, 4]

### TestAtomSelection.test_display_zero_atoms
_With no atoms selected the Add button is disabled and label shows 'None'._

- assert not dlg.add_button.isEnabled()
- assert 'None' in dlg.selection_label.text()

### TestAtomSelection.test_display_one_atom
_With one atom selected the Add button remains disabled._

- assert not dlg.add_button.isEnabled()

### TestAtomSelection.test_display_two_atoms_enables_add
_With two atoms selected the Add button is enabled and label shows 'Distance'._

- assert dlg.add_button.isEnabled()
- assert 'Distance' in dlg.selection_label.text()

### TestAtomSelection.test_display_three_atoms
_With three atoms selected the label shows 'Angle'._

- assert dlg.add_button.isEnabled()
- assert 'Angle' in dlg.selection_label.text()

### TestAtomSelection.test_display_four_atoms
_With four atoms selected the label shows 'Dihedral'._

- assert dlg.add_button.isEnabled()
- assert 'Dihedral' in dlg.selection_label.text()

### TestAddConstraint.test_distance_constraint_added
_Adding two atoms creates a Distance constraint with a positive value and default force._

- assert len(dlg.constraints) == 1
- assert ctype == 'Distance'
- assert cidx == (0, 1)
- assert cval > 0.0
- assert cforce == pytest.approx(100000.0)

### TestAddConstraint.test_distance_constraint_clears_selection
_add_constraint resets selected_atoms to an empty list after adding._

- assert dlg.selected_atoms == []

### TestAddConstraint.test_angle_constraint_added
_Adding three atoms creates an Angle constraint with value near tetrahedral angle._

- assert len(dlg.constraints) == 1
- assert ctype == 'Angle'
- assert cidx == (2, 0, 1)
- assert 90.0 <= cval <= 130.0

### TestAddConstraint.test_dihedral_constraint_added
_Adding four atoms creates a Dihedral constraint with correct atom indices._

- assert len(dlg.constraints) == 1
- assert ctype == 'Dihedral'
- assert cidx == (2, 0, 1, 5)

### TestAddConstraint.test_duplicate_constraint_rejected
_Adding the same constraint twice should not append a second entry._

- assert len(dlg.constraints) == 1

### TestAddConstraint.test_table_row_added_on_constraint
_add_constraint inserts a row into constraint_table._

- assert dlg.constraint_table.rowCount() == 1

### TestAddConstraint.test_custom_force_constant_used
_Force constant entered by the user is stored in the constraint tuple._

- assert cforce == pytest.approx(500.0)

### TestAddConstraint.test_invalid_force_constant_defaults
_Non-numeric force constant falls back to 1.0e5 (shows warning)._

- assert cforce == pytest.approx(100000.0)

### TestRemoveConstraint.test_remove_selected_row
_remove_constraint deletes the selected row and its corresponding constraint._

- assert len(dlg.constraints) == 1
- assert dlg.constraint_table.rowCount() == 1
- assert dlg.constraints[0][1] == (0, 2)

### TestRemoveConstraint.test_remove_all
_remove_all_constraints clears both the constraints list and the table._

- assert dlg.constraints == []
- assert dlg.constraint_table.rowCount() == 0

### TestRemoveConstraint.test_remove_all_on_empty_is_noop
_remove_all_constraints does not raise when called with no constraints._

- assert dlg.constraints == []

### TestRemoveConstraint.test_remove_button_disabled_after_remove_all
_Remove button is disabled after all constraints are cleared._

- assert not dlg.remove_button.isEnabled()

### TestCellChanged.test_edit_value_column_updates_constraint
_Editing the value cell updates the constraint's stored value._

- assert dlg.constraints[0][2] == pytest.approx(2.0)
- assert dlg.constraints[0][0] == original_type

### TestCellChanged.test_edit_force_column_updates_constraint
_Editing the force-constant cell updates the constraint's stored force value._

- assert dlg.constraints[0][3] == pytest.approx(25000.0)

### TestCellChanged.test_invalid_value_reverts_to_original
_Entering a non-numeric cell value reverts the cell to the original constraint value._

- assert dlg.constraints[0][2] == pytest.approx(original_val)

### TestCellChanged.test_non_value_column_ignored
_on_cell_changed does nothing when a non-editable column is touched._

- assert dlg.constraints[0] == original

### TestCellChanged.test_3element_constraint_upgraded_on_edit
_on_cell_changed must upgrade a legacy 3-element tuple to 4-element._

- assert len(dlg.constraints[0]) == 4
- assert dlg.constraints[0][2] == pytest.approx(1.8)
- assert dlg.constraints[0][3] == pytest.approx(100000.0)

### TestConstraintLoading.test_loads_4element_constraints
___init__ converts 4-element list constraints to tuples correctly._

- assert len(dlg.constraints) == 1
- assert ctype == 'Distance'
- assert cidx == (0, 1)
- assert cval == pytest.approx(1.54)
- assert cforce == pytest.approx(25000.0)

### TestConstraintLoading.test_loads_3element_constraints_with_default_force
_Legacy 3-element constraints should get default force 1.0e5._

- assert len(dlg.constraints) == 1
- assert ctype == 'Angle'
- assert cval == pytest.approx(109.5)
- assert cforce == pytest.approx(100000.0)

### TestConstraintLoading.test_loads_multiple_constraint_types
___init__ loads Distance, Angle, and Dihedral constraints simultaneously._

- assert len(dlg.constraints) == 3
- assert dlg.constraint_table.rowCount() == 3

### TestConstraintLoading.test_table_populated_on_load
_Preloaded constraints are shown in the table on dialog open._

- assert dlg.constraint_table.rowCount() == 1
- assert dlg.constraint_table.item(0, 0).text() == 'Distance'

### TestForceFieldMapping.test_ff_combo_set_from_opt_method
_Force-field combo is pre-selected based on the current optimization method._

- assert dlg.ff_combo.currentText() == expected_text

### TestForceFieldMapping.test_unknown_method_defaults_to_mmff94s
_An unrecognised optimization method defaults the force-field combo to MMFF94s._

- assert dlg.ff_combo.currentText() == 'MMFF94s'

### TestReject.test_reject_converts_tuples_to_lists
_reject serialises constraints as JSON-safe lists, not tuples._

- assert isinstance(saved, list)
- assert isinstance(saved[0], list)
- assert isinstance(saved[0][1], list)

### TestReject.test_reject_3element_constraint_gets_default_force
_Legacy 3-element tuples in constraints list must be serialised with default force._

- assert saved[0][3] == pytest.approx(100000.0)

### TestReject.test_reject_no_change_skips_state_update
_If saved constraints haven't changed, has_unsaved_changes must NOT be set._

- dlg.main_window.state_manager.update_window_title.assert_not_called()

### TestOptimizationExecution.test_apply_optimization_starts_thread
_apply_optimization disables the button and starts the worker thread._

- assert not dlg.optimize_button.isEnabled()
- mock_thread_class.assert_called_once()
- mock_thread.start.assert_called_once()

### TestOptimizationExecution.test_on_optimization_finished
__on_optimization_finished re-enables the button and redraws the molecule._

- assert dlg.optimize_button.isEnabled()
- dlg.main_window.view_3d_manager.draw_molecule_3d.assert_called_once_with(dlg.mol)
- dlg.main_window.edit_actions_manager.push_undo_state.assert_called_once()
- assert dlg.main_window.view_3d_manager.atom_positions_3d[0] == [1.0, 2.0, 3.0]

### TestOptimizationExecution.test_on_optimization_error
__on_optimization_error re-enables the button and shows a critical dialog._

- assert dlg.optimize_button.isEnabled()
- mock_mb.critical.assert_called_once()
- assert 'Test Error' in mock_mb.critical.call_args[0][2]

## tests/unit/test_custom_interactor_style.py

### test_custom_interactor_style_left_click_atom_selection
_Verify that left-clicking an atom in 3D edit mode successfully selects it._

- assert interactor_style._is_dragging_atom is True
- mock_parser_host.plotter.setCursor.assert_called_with(Qt.CursorShape.ClosedHandCursor)
- assert mock_parser_host.dragged_atom_info['id'] == 0

### test_custom_interactor_style_background_click
_Verify that clicking the background (no atom selected) allows VTK trackball rotation._

- assert interactor_style._is_dragging_atom is False
- mock_super_down.assert_called_once()

### test_measurement_atom_click_suppresses_vtk_release
_Atom measurement clicks consume press and must not leak release to VTK._

- mock_super_up.assert_not_called()

### test_atom_click_without_drag_resets_without_render_updates
_Simple atom clicks should not run drag redraw/update work._

- interactor_style.StopState.assert_called()
- mock_parser_host.edit_3d_manager.update_3d_selection_display.assert_not_called()
- mock_parser_host.edit_3d_manager.update_measurement_labels_display.assert_not_called()
- mock_parser_host.edit_3d_manager.update_2d_measurement_labels.assert_not_called()
- mock_parser_host.view_3d_manager.show_all_atom_info.assert_not_called()

### test_custom_interactor_style_move_selected_atoms_no_bfs
_Verify that clicking outside selection in MoveSelectedAtomsDialog toggles only that atom (no BFS)._

- mock_timer.assert_called_once()
- mock_dialog.on_atom_picked.assert_called_once_with(1)

### test_custom_interactor_style_right_click_rotation
_Verify right-click triggers group rotation and release finalizes it._

- assert mock_dialog.is_rotating_group_vtk is True
- assert mock_dialog.rotation_start_pos == (100, 100)
- assert mock_dialog.rotation_atom_idx == 0
- assert mock_dialog.group_centroid is not None
- assert mock_dialog.is_rotating_group_vtk is False
- assert mock_dialog.rotation_start_pos is None
- assert mock_dialog.is_rotating_group_vtk is False
- assert mock_dialog.rotation_start_pos is None
- mock_parser_host.view_3d_manager.draw_molecule_3d.assert_called_once()

## tests/unit/test_dialog_3d_picking_mixin.py

### test_init_defaults
_Dialog3DPickingMixin initialises with picking disabled and empty state._

- assert dlg.picking_enabled is False
- assert dlg._mouse_press_pos is None
- assert dlg._mouse_moved is False
- assert dlg.selection_labels == []

### test_eventfilter_none_event_returns_false
_eventFilter returns False immediately when the event is None._

- assert dlg.eventFilter(MagicMock(), None) is False

### test_eventfilter_plotter_none_returns_false
_eventFilter returns False when the plotter is None._

- assert dlg.eventFilter(MagicMock(), _left_press_event()) is False

### test_eventfilter_mol_none_returns_false
_eventFilter returns False when the molecule is None._

- assert dlg.eventFilter(MagicMock(), _left_press_event()) is False

### test_eventfilter_non_interactor_returns_false
_eventFilter returns False when the watched object is not the VTK interactor._

- assert dlg.eventFilter(other_obj, _left_press_event()) is False

### test_eventfilter_atom_click_calls_on_atom_picked
_A left-click on an atom calls on_atom_picked and returns True._

- assert result is True
- assert 0 in dlg._picked

### test_eventfilter_atom_click_consumes_matching_release
_The MouseButtonRelease following an atom pick is consumed (returns True)._

- assert dlg._mouse_press_pos is None
- assert dlg._mouse_moved is False
- assert dlg.eventFilter(plotter.interactor, _left_press_event()) is True
- assert dlg.eventFilter(plotter.interactor, _release_event()) is True

### test_eventfilter_atom_click_miss_returns_false
_A left-click that does not hit an atom actor returns False._

- assert result is False
- assert dlg._picked == []

### test_eventfilter_mouse_move_sets_moved_flag
_A mouse-move event beyond the drag threshold sets _mouse_moved to True._

- assert dlg._mouse_moved is True

### test_eventfilter_mouse_move_small_does_not_set_flag
_A mouse-move below the drag threshold does not set _mouse_moved._

- assert dlg._mouse_moved is False

### test_eventfilter_release_pure_click_calls_clear_selection
_A mouse-release with no drag calls clear_selection and resets press state._

- dlg.clear_selection.assert_called_once()
- assert dlg._mouse_press_pos is None

### test_eventfilter_release_after_drag_no_clear_selection
_A mouse-release after dragging does not call clear_selection._

- dlg.clear_selection.assert_not_called()
- assert dlg._mouse_press_pos is None

### test_enable_picking_installs_event_filter
_enable_picking installs the dialog as an event filter on the VTK interactor._

- plotter.interactor.installEventFilter.assert_called_once_with(dlg)
- assert dlg.picking_enabled is True

### test_enable_picking_none_plotter_no_crash
_enable_picking does not raise and leaves picking_enabled False when plotter is None._

- assert dlg.picking_enabled is False

### test_disable_picking_removes_event_filter
_disable_picking removes the event filter from the VTK interactor._

- plotter.interactor.removeEventFilter.assert_called_once_with(dlg)
- assert dlg.picking_enabled is False

### test_disable_picking_when_not_enabled_is_noop
_disable_picking does nothing when picking was not already enabled._

- plotter.interactor.removeEventFilter.assert_not_called()

### test_clear_atom_labels_removes_actors
_clear_atom_labels removes all label actors from the plotter and empties the list._

- plotter.remove_actor.assert_any_call(actor1)
- plotter.remove_actor.assert_any_call(actor2)
- assert dlg.selection_labels == []

### test_clear_atom_labels_none_plotter_empties_list
_clear_atom_labels empties selection_labels even when plotter is None._

- assert dlg.selection_labels == []

### test_add_selection_label_calls_add_point_labels
_add_selection_label calls plotter.add_point_labels and appends the actor._

- plotter.add_point_labels.assert_called_once()
- assert len(dlg.selection_labels) == 1

### test_add_selection_label_none_plotter_no_crash
_add_selection_label does not raise when plotter is None._

- assert dlg.selection_labels == []

### test_add_selection_label_none_positions_no_crash
_add_selection_label does not raise when atom_positions_3d is None._

- assert dlg.selection_labels == []

### test_show_atom_labels_for_clears_then_adds
_show_atom_labels_for clears existing labels then adds one per atom._

- plotter.remove_actor.assert_called()
- assert len(dlg.selection_labels) == 2

### test_show_atom_labels_for_empty_list_clears_all
_show_atom_labels_for with an empty list removes all existing labels._

- assert dlg.selection_labels == []

## tests/unit/test_dialog_logic.py

### test_bond_length_adjustment_logic
_Test the geometric logic of bond length adjustment directly._

- assert pytest.approx(final_dist, abs=0.001) == 2.0

### test_alignment_logic
_Test the geometry logic for aligning a bond to a specific axis._

- assert np.allclose(p0, [0, 0, 0], atol=1e-07)
- assert p1[0] > 0
- assert pytest.approx(p1[1], abs=1e-07) == 0
- assert pytest.approx(p1[2], abs=1e-07) == 0

### test_angle_adjustment_logic
_Test the geometric logic of bond angle adjustment directly._

- assert initial_angle != pytest.approx(120.0)
- assert pytest.approx(final_angle, abs=0.01) == 120.0

### test_dihedral_adjustment_logic
_Test the geometric logic of dihedral angle adjustment._

- assert pytest.approx(abs(final_dihedral), abs=0.01) == 180.0

### test_translation_logic
_Test the geometric logic of centroid-based translation._

- assert np.allclose(new_pos0, expected_pos0, atol=1e-07)
- assert np.allclose(new_pos5, p5_initial, atol=1e-07)

### test_move_group_logic
_Test the translation and rotation logic in MoveGroupDialog._

- assert np.allclose(np.array(mol.GetConformer().GetAtomPosition(0)), initial_pos0 + [5, 0, 0])
- assert np.allclose(np.array(mol.GetConformer().GetAtomPosition(5)), initial_pos5)
- assert np.allclose(new_pos0, expected_rotated, atol=1e-07)

### test_abs_multi_atom_selection
_Clicking a second atom in Absolute tab accumulates the selection._

- assert dialog.selected_atoms == {0}
- assert dialog.selected_atoms == {0, 3}

### test_abs_coordinate_inputs_populated_on_pick
_Selecting an atom auto-fills X/Y/Z inputs with its current position._

- dialog.abs_x_input.setText.assert_called_once_with(f'{pos.x:.4f}')
- dialog.abs_y_input.setText.assert_called_once_with(f'{pos.y:.4f}')
- dialog.abs_z_input.setText.assert_called_once_with(f'{pos.z:.4f}')

### test_abs_apply_moves_entire_molecule_by_default
_apply_absolute with move_mol=True shifts every atom by the same delta._

- assert np.allclose(new_pos0, atom0_pos + [5, 0, 0], atol=0.0001)
- assert np.allclose(new_pos5, atom5_pos + [5, 0, 0], atol=0.0001)

### test_abs_apply_moves_only_selected_atom
_apply_absolute with move_mol=False moves only the selected atom._

- assert np.allclose(new_pos0, atom0_pos + [0, 3, 0], atol=0.0001)
- assert np.allclose(new_pos5, atom5_pos, atol=0.0001)

### test_abs_apply_noop_when_no_delta
_apply_absolute does nothing when target equals current position._

- assert np.allclose(mol.GetConformer().GetPositions(), original_positions, atol=1e-07)

### test_abs_apply_requires_one_atom
_apply_absolute shows a warning and does nothing when no atom is selected._

- assert np.allclose(mol.GetConformer().GetPositions(), original_positions, atol=1e-07)
- mock_mb.warning.assert_called_once()

### test_tab_switch_clears_selection
_Switching tabs clears selected_atoms and atom labels._

- assert dialog.selected_atoms == set()
- mock_clear.assert_called_once()

### test_delta_tab_atom_pick_toggles
_In Delta tab, picking an already-selected atom deselects it._

- assert 0 not in dialog.selected_atoms

## tests/unit/test_dialog_manager.py

### TestGetPreselectedAtoms3D.test_returns_empty_when_no_selection
__get_preselected_atoms_3d returns [] when no atoms are selected._

- assert dm._get_preselected_atoms_3d() == []

### TestGetPreselectedAtoms3D.test_returns_selected_atoms
__get_preselected_atoms_3d returns the current measurement selection._

- assert result == [1, 2, 3]

### TestGetPreselectedAtoms3D.test_logs_error_when_edit_3d_manager_missing
__get_preselected_atoms_3d returns [] and logs when edit_3d_manager is absent._

- assert result == []

### TestShowAboutDialog.test_creates_and_execs_dialog
_show_about_dialog creates an AboutDialog and calls exec._

- MockAbout.assert_called_once_with(dm.host, dm.host)
- instance.exec.assert_called_once()

### TestOpenPeriodicTableDialog.test_creates_connects_and_execs
_open_periodic_table_dialog creates, connects element_selected, and execs the dialog._

- MockPT.assert_called_once_with(dm.host)
- instance.element_selected.connect.assert_called_once_with(dm.host.ui_manager.set_atom_from_periodic_table)
- instance.exec.assert_called_once()

### TestOpenPeriodicTableDialog.test_unchecks_tool_group_action
_open_periodic_table_dialog unchecks the active tool-group action._

- checked.setChecked.assert_called_with(False)

### TestOpenAnalysisWindow.test_opens_when_mol_exists
_open_analysis_window creates and execs AnalysisWindow when a molecule is loaded._

- MockAW.assert_called_once_with(dm.host.view_3d_manager.current_mol, dm.host, is_xyz_derived=dm.host.is_xyz_derived)
- instance.exec.assert_called_once()

### TestOpenAnalysisWindow.test_shows_error_when_no_mol
_open_analysis_window shows a status error when no 3D molecule is loaded._

- MockAW.assert_not_called()
- dm.host.statusBar_mock.showMessage.assert_called_once()
- assert '3D' in msg or 'generate' in msg.lower()

### TestOpenTemplateDialog.test_creates_and_execs
_open_template_dialog creates a UserTemplateDialog and calls exec._

- MockUT.assert_called_once_with(dm.host, dm.host)
- instance.exec.assert_called_once()

### TestOpenTemplateDialogAndActivate.test_creates_new_dialog_when_none_exists
_open_template_dialog_and_activate creates a new modeless dialog when none is open._

- MockUT.assert_called_once_with(dm.host, dm.host)
- instance.show.assert_called_once()
- instance.finished.connect.assert_called_once()

### TestOpenTemplateDialogAndActivate.test_raises_existing_visible_dialog
_open_template_dialog_and_activate raises and activates an already-open dialog._

- MockUT.assert_not_called()
- existing.raise_.assert_called_once()
- existing.activateWindow.assert_called_once()

### TestOpenTemplateDialogAndActivate.test_on_finished_sets_mode_when_template_selected
_Finishing with a template selected activates template mode and shows a status message._

- assert captured_cb
- dm.host.ui_manager.set_mode.assert_called_once_with('template_user_benzene')
- dm.host.statusBar_mock.showMessage.assert_called_once()

### TestOpenTemplateDialogAndActivate.test_on_finished_noop_when_no_template_selected
_Finishing without a selection does not change mode._

- dm.host.ui_manager.set_mode.assert_not_called()

### TestSave2DAsTemplate.test_warns_when_no_atoms
_save_2d_as_template shows a warning when the canvas has no atoms._

- mock_warn.assert_called_once()
- assert 'No structure' in args[2] or 'template' in args[2].lower()

### TestSave2DAsTemplate.test_noop_on_cancelled_input
_save_2d_as_template does nothing when the name dialog is cancelled._

- dm.host.state_manager.data.to_template_dict.assert_not_called()

### TestSave2DAsTemplate.test_noop_on_blank_name
_save_2d_as_template does nothing when a blank name is entered._

- dm.host.state_manager.data.to_template_dict.assert_not_called()

### TestSave2DAsTemplate.test_saves_template_file
_save_2d_as_template writes a .pmetmplt file to the user-templates directory._

- assert saved.exists()
- assert json.loads(saved.read_text())['name'] == 'mytemplate'

### TestSave2DAsTemplate.test_overwrites_after_yes_confirmation
_save_2d_as_template overwrites an existing file when the user confirms Yes._

- assert 'atoms' in json.loads(f.read_text())

### TestSave2DAsTemplate.test_skips_overwrite_on_no_confirmation
_save_2d_as_template leaves the existing file unchanged when user answers No._

- assert json.loads(f.read_text()) == {'original': True}

### TestSave2DAsTemplate.test_shows_error_on_exception
_save_2d_as_template shows a critical error dialog when serialisation raises._

- mock_crit.assert_called_once()

### TestModelessGeometryDialogs.test_open_translation_dialog
_open_translation_dialog shows a modeless dialog and wires up all signals._

- _assert_modeless(dm, 'open_translation_dialog', 'TranslationDialog')

### TestModelessGeometryDialogs.test_translation_disables_measurement_mode
_open_translation_dialog disables measurement mode before showing the dialog._

- host.edit_3d_manager.toggle_measurement_mode.assert_called_with(False)

### TestModelessGeometryDialogs.test_open_move_group_dialog
_open_move_group_dialog shows a modeless dialog and wires up all signals._

- _assert_modeless(dm, 'open_move_group_dialog', 'MoveGroupDialog')

### TestModelessGeometryDialogs.test_open_align_plane_dialog
_open_align_plane_dialog shows a modeless dialog and wires up all signals._

- _assert_modeless(dm, 'open_align_plane_dialog', 'AlignPlaneDialog', 'xy')

### TestModelessGeometryDialogs.test_align_plane_message_contains_plane
_The accepted callback for align-plane includes the plane name in the status message._

- assert 'XZ' in host.statusBar_mock.showMessage.call_args[0][0]

### TestModelessGeometryDialogs.test_open_planarize_dialog
_open_planarize_dialog shows a modeless dialog and wires up all signals._

- _assert_modeless(dm, 'open_planarize_dialog', 'PlanarizeDialog')

### TestModelessGeometryDialogs.test_open_alignment_dialog
_open_alignment_dialog shows a modeless dialog and wires up all signals._

- _assert_modeless(dm, 'open_alignment_dialog', 'AlignmentDialog', 'x')

### TestModelessGeometryDialogs.test_alignment_message_contains_axis
_The accepted callback for alignment includes the axis name in the status message._

- assert 'Y' in host.statusBar_mock.showMessage.call_args[0][0]

### TestModelessGeometryDialogs.test_open_bond_length_dialog
_open_bond_length_dialog shows a modeless dialog and wires up all signals._

- _assert_modeless(dm, 'open_bond_length_dialog', 'BondLengthDialog')

### TestModelessGeometryDialogs.test_open_angle_dialog
_open_angle_dialog shows a modeless dialog and wires up all signals._

- _assert_modeless(dm, 'open_angle_dialog', 'AngleDialog')

### TestModelessGeometryDialogs.test_open_dihedral_dialog
_open_dihedral_dialog shows a modeless dialog and wires up all signals._

- _assert_modeless(dm, 'open_dihedral_dialog', 'DihedralDialog')

### TestModelessGeometryDialogs.test_accepted_status_messages
_Each dialog's first accepted lambda posts the right status bar message._

- assert actual == expected

### TestModelessGeometryDialogs.test_accepted_pushes_undo_state
_Second accepted lambda calls push_undo_state._

- host.edit_actions_manager.push_undo_state.assert_called_once()

### TestModelessGeometryDialogs.test_finished_removes_dialog_from_list
_finished lambda calls remove_dialog_from_list with this dialog._

- host.edit_3d_manager.remove_dialog_from_list.assert_called_once_with(instance)

### TestOpenMirrorDialog.test_opens_when_mol_exists
_open_mirror_dialog shows a modeless MirrorDialog when a molecule is loaded._

- MockM.assert_called_once_with(dm.host.view_3d_manager.current_mol, dm.host, parent=dm.host)
- instance.show.assert_called_once()
- instance.finished.connect.assert_called_once()
- assert instance in dm.host.edit_3d_manager.active_3d_dialogs

### TestOpenMirrorDialog.test_shows_error_when_no_mol
_open_mirror_dialog shows a status error when no 3D molecule is loaded._

- MockM.assert_not_called()
- dm.host.statusBar_mock.showMessage.assert_called_with('No 3D molecule loaded.')

### TestOpenMirrorDialog.test_disables_measurement_mode
_open_mirror_dialog disables measurement mode before showing the dialog._

- host.edit_3d_manager.toggle_measurement_mode.assert_called_with(False)

### TestOpenSettingsDialog.test_creates_and_execs
_open_settings_dialog creates a SettingsDialog and calls exec._

- MockSD.assert_called_once_with(dm.host.init_manager.settings, parent=dm.host)
- instance.exec.assert_called_once()

### TestOpenColorSettingsDialog.test_creates_and_execs
_open_color_settings_dialog creates a ColorSettingsDialog and calls exec._

- MockCD.assert_called_once_with(dm.host.init_manager.settings, parent=dm.host)
- instance.exec.assert_called_once()

### TestOpenConstrainedOptimizationDialog.test_opens_when_mol_exists
_open_constrained_optimization_dialog shows a modeless dialog when a molecule is loaded._

- MockCO.assert_called_once_with(dm.host.view_3d_manager.current_mol, dm.host, parent=dm.host)
- instance.show.assert_called_once()
- instance.finished.connect.assert_called_once()
- assert instance in dm.host.edit_3d_manager.active_3d_dialogs

### TestOpenConstrainedOptimizationDialog.test_shows_error_when_no_mol
_open_constrained_optimization_dialog shows a status error when no molecule is loaded._

- MockCO.assert_not_called()
- dm.host.statusBar_mock.showMessage.assert_called_with('No 3D molecule loaded.')

### TestOpenConstrainedOptimizationDialog.test_disables_measurement_mode
_open_constrained_optimization_dialog disables measurement mode before showing._

- host.edit_3d_manager.toggle_measurement_mode.assert_called_with(False)

### TestOpenConstrainedOptimizationDialog.test_finished_removes_from_active_dialogs
_finished callback removes the dialog from the active dialogs list._

- dm.host.edit_3d_manager.remove_dialog_from_list.assert_called_once_with(instance)

## tests/unit/test_edit_3d_logic.py

### test_calculate_distance_logic
_Verify calculation of distance between 3D points._

- assert dist == pytest.approx(1.5)

### test_calculate_angle_logic
_Verify calculation of angle between three 3D points._

- assert angle == pytest.approx(90.0)

### test_calculate_dihedral_logic
_Verify calculation of dihedral angle between four 3D points._

- assert abs(dihedral - 45.0) < 0.01

### test_handle_measurement_atom_selection
_Verify handling of atom selection for measurements._

- assert 10 in edit3d.selected_atoms_for_measurement
- assert len(edit3d.selected_atoms_for_measurement) == 0

### test_clear_3d_selection
_Verify clearing of 3D selection highlights._

- assert edit3d.plotter.remove_actor.called

### test_toggle_measurement_mode
_Verify toggling of measurement mode._

- assert edit3d.measurement_mode is True
- assert edit3d.statusBar().showMessage.called
- assert edit3d.measurement_mode is False

### test_calculate_and_display_measurements_trigger
_Verify triggering of measurement calculation and display._

- assert mock_display.called
- assert 'Distance' in mock_display.call_args[0][0][0]

### test_add_2d_measurement_label_adds_item
_Verify add_2d_measurement_label creates a label and registers it._

- assert len(edit3d.measurement_label_items_2d) == 1
- assert isinstance(label, QGraphicsTextItem)
- assert label.toPlainText() == '1.234 Å'
- assert label.zValue() == 2000
- mock_parser_host.init_manager.scene.addItem.assert_called_once_with(label)
- assert pos.x() == pytest.approx(16.0)
- assert pos.y() == pytest.approx(8.0)

### test_add_2d_measurement_label_accumulates
_Verify multiple labels are each appended to measurement_label_items_2d._

- assert len(edit3d.measurement_label_items_2d) == 3
- assert mock_parser_host.init_manager.scene.addItem.call_count == 3

### test_clear_2d_measurement_labels_removes_items
_Verify clear_2d_measurement_labels removes all items from the scene._

- mock_parser_host.init_manager.scene.removeItem.assert_called_once_with(mock_label)
- assert edit3d.measurement_label_items_2d == []

### test_clear_2d_measurement_labels_skips_deleted
_Verify clear_2d_measurement_labels skips sip-deleted items._

- mock_parser_host.init_manager.scene.removeItem.assert_not_called()
- assert edit3d.measurement_label_items_2d == []

### test_clear_2d_measurement_labels_skips_unscened
_Verify that items not attached to a scene are not passed to removeItem._

- mock_parser_host.init_manager.scene.removeItem.assert_not_called()
- assert edit3d.measurement_label_items_2d == []

### test_find_rdkit_atom_index_no_mol
_Verify find_rdkit_atom_index returns None when no molecule is loaded._

- assert edit3d.find_rdkit_atom_index(atom_item) is None

### test_find_rdkit_atom_index_no_item
_Verify find_rdkit_atom_index returns None for a None atom_item._

- assert edit3d.find_rdkit_atom_index(None) is None

### test_find_rdkit_atom_index_with_map
_Verify find_rdkit_atom_index returns the mapped RDKit index._

- assert edit3d.find_rdkit_atom_index(atom_item) == 2

### test_find_rdkit_atom_index_missing_map_entry
_Verify find_rdkit_atom_index returns None when atom_id not in map._

- assert edit3d.find_rdkit_atom_index(atom_item) is None

### test_display_measurement_text_empty_clears_actor
_Verify display_measurement_text with empty list removes existing actor._

- mock_parser_host.view_3d_manager.plotter.remove_actor.assert_called_once()
- assert edit3d.measurement_text_actor is None

### test_display_measurement_text_adds_actor
_Verify display_measurement_text calls add_text with joined lines._

- assert 'Distance 1-2: 1.540 Å' in text_arg
- assert 'Angle: 109.5°' in text_arg
- assert call_kwargs[1].get('color') == 'white'

### test_display_measurement_text_light_background
_Verify display_measurement_text uses black text on a light background._

- assert call_kwargs[1].get('color') == 'black'

### test_toggle_atom_selection_3d_add
_Verify toggle_atom_selection_3d adds an atom not yet selected._

- assert 3 in edit3d.selected_atoms_3d
- mock_update.assert_called_once()

### test_toggle_atom_selection_3d_remove
_Verify toggle_atom_selection_3d removes an already-selected atom._

- assert 3 not in edit3d.selected_atoms_3d
- mock_update.assert_called_once()

### test_toggle_atom_selection_3d_idempotent_add
_Verify toggling different atoms accumulates them independently._

- assert edit3d.selected_atoms_3d == {1, 2}

### test_toggle_on_when_edit_mode_active_disables_edit_mode
_Enabling measurement mode while 3D edit mode is active disables edit mode._

- host.init_manager.edit_3d_action.setChecked.assert_called_once_with(False)
- host.ui_manager.toggle_3d_edit_mode.assert_called_once_with(False)

### test_toggle_on_closes_active_dialogs
_Enabling measurement mode closes any open 3D edit dialogs._

- dlg.close.assert_called_once()
- assert mgr.active_3d_dialogs == []

### test_close_all_closes_each_dialog
_close_all_3d_edit_dialogs closes every dialog and empties the list._

- d1.close.assert_called_once()
- d2.close.assert_called_once()
- assert mgr.active_3d_dialogs == []

### test_close_all_handles_close_error
_close_all_3d_edit_dialogs swallows exceptions from dialog.close()._

- assert mgr.active_3d_dialogs == []

### test_close_all_empty_list_is_noop
_close_all_3d_edit_dialogs with an empty list does not raise._


### test_update_labels_display_adds_point_labels
_update_measurement_labels_display calls add_point_labels when labels exist._

- host.view_3d_manager.plotter.add_point_labels.assert_called_once()

### test_update_labels_display_no_labels_returns_early
_update_measurement_labels_display skips add_point_labels when list is empty._

- host.view_3d_manager.plotter.add_point_labels.assert_not_called()

### test_update_labels_display_no_mol_returns_early
_update_measurement_labels_display skips rendering when no molecule is loaded._

- host.view_3d_manager.plotter.add_point_labels.assert_not_called()

### test_clear_measurement_selection_clears_state
_clear_measurement_selection empties selected atoms and measurement labels._

- assert mgr.selected_atoms_for_measurement == []
- assert mgr.measurement_labels == []
- host.view_3d_manager.plotter.render.assert_called()

### test_clear_measurement_selection_removes_text_actor
_clear_measurement_selection removes the measurement text actor from the plotter._

- host.view_3d_manager.plotter.remove_actor.assert_any_call(actor)
- assert mgr.measurement_text_actor is None

### test_update_2d_labels_no_mol_returns_early
_update_2d_measurement_labels is a no-op when no molecule is loaded._

- mock_add.assert_not_called()

### test_update_2d_labels_no_atoms_data_returns_early
_update_2d_measurement_labels skips rendering when atoms data dict is empty._

- mock_add.assert_not_called()

### test_update_2d_labels_maps_atom_and_adds_label
_update_2d_measurement_labels resolves atom items and calls add_2d_measurement_label._

- mock_add.assert_called_once_with(atom_item, '1')

### test_update_3d_selection_empty_renders
_update_3d_selection_display renders even with an empty selection set._

- host.view_3d_manager.plotter.render.assert_called()

### test_update_3d_selection_no_mol_renders
_update_3d_selection_display still calls render when no molecule is loaded._

- host.view_3d_manager.plotter.render.assert_called()

### test_remove_dialog_present
_remove_dialog_from_list removes the dialog when it is in active_3d_dialogs._

- assert dlg not in mgr.active_3d_dialogs

### test_remove_dialog_absent_is_noop
_remove_dialog_from_list is a no-op when the dialog is not in the list._


### test_calculate_and_display_3_atoms_includes_angle
_Three selected atoms produces both a distance and an angle measurement._

- assert any(('Angle' in l for l in lines))
- assert any(('Distance' in l for l in lines))

### test_calculate_and_display_4_atoms_includes_dihedral
_Four selected atoms produces distance, angle, and dihedral measurements._

- assert any(('Dihedral' in l for l in lines))
- assert any(('Angle' in l for l in lines))
- assert any(('Distance' in l for l in lines))

### test_calculate_and_display_1_atom_does_nothing
_One selected atom does not trigger any measurement display._

- mock_disp.assert_not_called()

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

### TestEditActionsExtended.test_apply_chem_check_force_skip
_force_skip=True bypasses the chemistry check without setting any flags._

- assert manager.host.chem_check_tried is False
- assert manager.host.chem_check_failed is False

### TestEditActionsExtended.test_apply_chem_check_settings_skip
_skip_chemistry_checks=True in settings prevents the check from running._

- assert manager.host.chem_check_tried is False
- assert manager.host.chem_check_failed is False

### TestEditActionsExtended.test_apply_chem_check_success
_A valid molecule sets chem_check_tried=True and chem_check_failed=False._

- assert manager.host.chem_check_tried is True
- assert manager.host.chem_check_failed is False

### TestEditActionsExtended.test_apply_chem_check_failure
_A sanitization error sets chem_check_failed=True and disables the 3D optimize button._

- assert manager.host.chem_check_tried is True
- assert manager.host.chem_check_failed is True
- assert manager.host.init_manager.optimize_3d_button.setEnabled.called
- manager.host.init_manager.optimize_3d_button.setEnabled.assert_called_with(False)

### TestEditActionsExtended.test_clear_xyz_flags_with_mol_arg
__clear_xyz_flags removes XYZ-specific props and clears is_xyz_derived._

- assert not mol.HasProp('_xyz_skip_checks')
- assert not hasattr(mol, '_xyz_skip_checks')
- assert not hasattr(mol, 'xyz_atom_data')
- assert manager.host.is_xyz_derived is False
- manager.host.init_manager.optimize_3d_button.setEnabled.assert_called_with(True)

### TestEditActionsExtended.test_update_edit_menu_actions
_update_edit_menu_actions enables cut/copy/paste when items are selected._

- manager.host.init_manager.cut_action.setEnabled.assert_called_with(True)
- manager.host.init_manager.copy_action.setEnabled.assert_called_with(True)
- manager.host.init_manager.paste_action.setEnabled.assert_called_with(True)

### TestEditActionsExtended.test_open_rotate_2d_dialog
_open_rotate_2d_dialog calls rotate_molecule_2d with the user's chosen angle._

- manager.rotate_molecule_2d.assert_called_with(45.0)
- assert manager.last_rotation_angle == 45.0

### TestEditActionsExtended.test_rotate_molecule_2d_full
_rotate_molecule_2d updates atom positions and stores them in molecule data._

- atom1.setPos.assert_called()
- manager.host.state_manager.data.set_atom_pos.assert_called()

### TestEditActionsExtended.test_select_all
_select_all selects every AtomItem and BondItem in the scene._

- atom.setSelected.assert_called_with(True)
- bond.setSelected.assert_called_with(True)

### TestEditActionsExtended.test_clear_all
_clear_all(skip_check=True) returns True and shows the cleared status message._

- assert result is True
- manager.host.statusBar().showMessage.assert_called_with('Cleared all data.')

### TestEditActionsExtended.test_cut_selection
_cut_selection copies then deletes the selected items._

- manager.copy_selection.assert_called()
- manager.host.init_manager.scene.delete_items.assert_called()
- manager.host.statusBar().showMessage.assert_called_with('Cut selection.', 2000)

### TestEditActionsExtended.test_cut_selection_no_selection
_cut_selection is a no-op when nothing is selected._

- manager.copy_selection.assert_not_called()

### TestEditActionsExtended.test_adjust_molecule_positions_no_collision
_Fragments already far apart are not moved by collision avoidance._

- assert list(conf.GetAtomPosition(0)) == pos0_before
- assert list(conf.GetAtomPosition(1)) == pos1_before

### TestEditActionsExtended.test_adjust_molecule_positions_with_collision
_Colliding fragments are separated by collision avoidance._

- assert not np.array_equal(pos0_before, pos0_after)
- assert not np.array_equal(pos1_before, pos1_after)
- assert np.linalg.norm(pos0_after - pos1_after) > np.linalg.norm(pos0_before - pos1_before)

### TestEditActionsExtended.test_adjust_molecule_positions_single_fragment
_A single fragment is not moved by collision avoidance._

- assert list(conf.GetAtomPosition(0)) == pos_before

### TestEditActionsExtended.test_apply_chem_check_missing_button
_apply_chem_check_and_set_flags does not raise when optimize_3d_button is absent._


### TestEditActionsExtended.test_clear_xyz_flags_current_mol
__clear_xyz_flags(mol=None) uses the current_mol from the view manager._

- assert not mol.HasProp('_xyz_skip_checks')
- assert manager.host.is_xyz_derived is False

### TestEditActionsExtended.test_clear_xyz_flags_missing_zoom
__clear_xyz_flags does not raise when reset_zoom is absent._


## tests/unit/test_export_logic.py

### test_create_multi_material_obj_advanced
_Verify that create_multi_material_obj creates both .obj and .mtl files._

- assert os.path.exists(obj_path)
- assert os.path.exists(mtl_path)
- assert os.path.getsize(obj_path) > 0
- assert os.path.getsize(mtl_path) > 0
- assert 'mtllib' in content
- assert 'v ' in content

### test_export_2d_png_basic_trigger
_Verify that export_2d_png successfully creates a PNG file after user confirmation._

- assert os.path.exists(save_path)
- assert os.path.getsize(save_path) > 0

### test_export_2d_svg_trigger
_Verify that export_2d_svg successfully creates an SVG file._

- assert os.path.exists(save_path)
- assert os.path.getsize(save_path) > 0
- assert '<svg' in f.read().lower()

### test_export_stl_allows_missing_current_mol
_Verify that export_stl can proceed when plugins do not set current_mol._

- assert mock_file_dialog.called
- assert not any((old_error in str(args) for args, _ in exporter.statusBar().showMessage.call_args_list))

### test_export_obj_mtl_allows_missing_current_mol
_Verify that export_obj_mtl can proceed when plugins do not set current_mol._

- assert mock_file_dialog.called
- assert not any((old_error in str(args) for args, _ in exporter.statusBar().showMessage.call_args_list))

### test_export_stl_success_trigger
_Verify that export_stl triggers the mesh save logic for a valid 3D molecule._

- assert args[0] == save_path
- assert kwargs.get('binary') is True
- assert mesh.save.called

### test_export_obj_mtl_success_trigger
_Verify that export_obj_mtl triggers the multi-material OBJ creation logic._

- assert mock_file_dialog.called
- assert mock_create.called

### test_export_3d_png_logic
_Verify that export_3d_png triggers the plotter screenshot logic._

- assert mock_plotter.screenshot.called

### test_export_3d_png_allows_missing_current_mol
_Verify that export_3d_png can proceed when plugins do not set current_mol._

- assert mock_file_dialog.called
- assert not any((old_error in str(args) for args, _ in exporter.statusBar().showMessage.call_args_list))

### test_export_color_stl_logic
_Verify that export_color_stl triggers the status bar message (success indicator)._

- assert exporter.statusBar().showMessage.called

### test_export_color_stl_allows_missing_current_mol
_Verify that export_color_stl can proceed when plugins do not set current_mol._

- assert mock_file_dialog.called
- assert not any((old_error in str(args) for args, _ in exporter.statusBar().showMessage.call_args_list))

### test_export_from_3d_view_with_colors_complex_splitting
_Test the complex logic of splitting a mesh by per-vertex colors._

- assert len(res) == 2
- assert colors[0] == (0, 0, 255)
- assert colors[1] == (255, 0, 0)

### test_export_2d_png_hides_items
_Test that export_2d_png hides non-atom items and restores them._

- assert other_item.hide.called
- other_item.setVisible.assert_called_with(True)

## tests/unit/test_geometry.py

### test_3d_bond_lengths
_Verify optimized 3D coordinates yield physical bond lengths._

- assert mol is not None
- assert cc_bond is not None
- assert 1.51 < dist < 1.57
- assert 107 < angle < 112

### test_mirror_dialog_logic
_Verify MirrorDialog correctly manipulates coordinates and UI (software logic)._

- assert x_after == pytest.approx(-x_before, abs=0.001)
- assert main_window.view_3d_manager.draw_molecule_3d.called
- assert main_window.view_3d_manager.update_chiral_labels.called
- assert main_window.edit_actions_manager.push_undo_state.called

### test_planarize_logic
_Verify planarize functionality using the actual PlanarizeDialog logic._

- assert s[-1] < 1e-10
- assert main_window.view_3d_manager.draw_molecule_3d.called
- assert main_window.view_3d_manager.update_chiral_labels.called
- assert main_window.edit_actions_manager.push_undo_state.called

### test_rodrigues_rotate_90_deg
_Rotate [1,0,0] by 90° around [0,0,1] → expect [0,1,0]._

- np.testing.assert_allclose(result, [0.0, 1.0, 0.0], atol=1e-12)

### test_rodrigues_rotate_identity
_Rotate by 0° → vector unchanged._

- np.testing.assert_allclose(result, v, atol=1e-12)

### test_adjust_bond_angle_simple
_Set a 90° angle to 120° and verify._

- assert angle_deg == pytest.approx(target, abs=1e-08)
- np.testing.assert_allclose(positions[1], [0, 0, 0], atol=1e-12)
- assert np.linalg.norm(positions[2] - positions[1]) == pytest.approx(1.0, abs=1e-12)

### test_adjust_bond_angle_with_group
_Move multiple atoms; verify relative geometry is preserved._

- assert angle_deg == pytest.approx(60.0, abs=1e-08)
- assert cd_after == pytest.approx(cd_before, abs=1e-12)

### test_adjust_bond_angle_collinear
_Collinear atoms → fallback axis is used, rotation still succeeds._

- assert delta != 0.0
- assert angle_deg == pytest.approx(target, abs=1e-08)
- assert np.linalg.norm(positions[2] - positions[1]) == pytest.approx(1.0, abs=1e-12)

### test_optimize_2d_coords
_Verify 2D coordinate optimization generating coordinates for a simple molecule._

- assert len(new_pos) == 6
- assert len(pos) == 2

### test_calculate_best_fit_plane_projection
_Verify orthogonal projection onto a best-fit plane._

- assert ps[-1] < 1e-10

### test_rotate_2d_points
_Verify 2D rotation of point maps._

- np.testing.assert_allclose(rotated[1], [0, 1], atol=1e-12)
- np.testing.assert_allclose(rotated[2], [-1, 0], atol=1e-12)

### test_resolve_2d_overlaps
_Verify 2D overlap resolution logic handles collisions correctly._

- assert len(moves) > 0
- assert len(moves[0][0]) == 1

## tests/unit/test_geometry_base_dialog.py

### TestSyncInputToSlider.test_valid_float_sets_slider
_A valid float string sets the slider to the scaled integer value._

- assert s.value() == 250

### TestSyncInputToSlider.test_invalid_text_is_ignored
_Non-numeric input leaves the slider value unchanged._

- assert s.value() == 123

### TestSyncInputToSlider.test_wrap_true_normalises_into_range
_wrap=True normalises angles outside (-180, 180] into that range._

- assert s.value() == -90

### TestSyncInputToSlider.test_wrap_false_uses_raw_value
_wrap=False sets the slider to the raw float value without normalising._

- assert s.value() == 90

### TestOnSliderPressed.test_incomplete_selection_is_noop
_on_slider_pressed is a no-op when selection is incomplete._

- assert not dlg._slider_dragging
- assert dlg._snapshot_positions is None

### TestOnSliderPressed.test_complete_selection_sets_dragging_and_snapshot
_on_slider_pressed with complete selection sets _slider_dragging and snapshots positions._

- assert dlg._slider_dragging is True
- assert dlg._snapshot_positions is not None
- assert len(dlg._snapshot_positions) == dlg.mol.GetNumAtoms()

### TestOnSliderReleased.test_clears_dragging_flag
_on_slider_released sets _slider_dragging to False._

- assert dlg._slider_dragging is False

### TestOnSliderReleased.test_calls_draw_molecule_3d
_on_slider_released triggers a 3D redraw of the molecule._

- dlg.main_window.view_3d_manager.draw_molecule_3d.assert_called_once_with(dlg.mol)

### TestOnSliderValueChangedClick.test_skips_when_dragging
_on_slider_value_changed_click returns early while the slider is being dragged._

- assert dlg.applied_values == []

### TestOnSliderValueChangedClick.test_skips_when_selection_incomplete
_on_slider_value_changed_click returns early when selection is incomplete._

- assert dlg.applied_values == []

### TestOnSliderValueChangedClick.test_updates_input_and_calls_apply
_on_slider_value_changed_click updates the text box and calls apply_geometry_update._

- assert inp.text() == '1.540'
- assert dlg.applied_values == [pytest.approx(1.54)]

### TestOnSliderMovedRealtime.test_skips_when_incomplete
_on_slider_moved_realtime is a no-op when selection is incomplete._

- assert dlg.applied_values == []

### TestOnSliderMovedRealtime.test_updates_input_and_calls_apply
_on_slider_moved_realtime updates the text box and calls apply_geometry_update._

- assert inp.text() == '2.000'
- assert dlg.applied_values == [pytest.approx(2.0)]

### TestLoggingErrorFallback.test_missing_chiral_labels_on_released
_on_slider_released does not raise when update_chiral_labels is absent._


### TestLoggingErrorFallback.test_missing_chiral_labels_on_click
_on_slider_value_changed_click does not raise when update_chiral_labels is absent._


## tests/unit/test_geometry_dialogs.py

### TestBondLengthPicking.test_first_pick_sets_atom1
_First atom pick sets atom1_idx and leaves atom2_idx None._

- assert dlg.atom1_idx == 0
- assert dlg.atom2_idx is None

### TestBondLengthPicking.test_second_pick_sets_atom2
_Second atom pick fills atom2_idx._

- assert dlg.atom1_idx == 0
- assert dlg.atom2_idx == 1

### TestBondLengthPicking.test_third_pick_resets_to_new_atom1
_A third pick restarts selection with the new atom as atom1._

- assert dlg.atom1_idx == 2
- assert dlg.atom2_idx is None

### TestBondLengthPicking.test_clear_selection
_clear_selection resets both atom indices to None._

- assert dlg.atom1_idx is None
- assert dlg.atom2_idx is None

### TestBondLengthPicking.test_preselected_atoms_loaded
_BondLengthDialog pre-populates atom indices from the preselected_atoms list._

- assert dlg.atom1_idx == 0
- assert dlg.atom2_idx == 1

### TestBondLengthIsComplete.test_incomplete_with_one_atom
__is_selection_complete returns False with only one atom selected._

- assert not dlg._is_selection_complete()

### TestBondLengthIsComplete.test_complete_with_two_atoms
__is_selection_complete returns True when both atom indices are set._

- assert dlg._is_selection_complete()

### TestBondLengthUpdateDisplay.test_no_atoms_disables_apply
_update_display disables the Apply button and shows 'No atoms' with no selection._

- assert not dlg.apply_button.isEnabled()
- assert 'No atoms' in dlg.selection_label.text()

### TestBondLengthUpdateDisplay.test_one_atom_shows_symbol
_update_display shows the first atom's symbol and keeps Apply disabled._

- assert not dlg.apply_button.isEnabled()
- assert mol.GetAtomWithIdx(0).GetSymbol() in dlg.selection_label.text()

### TestBondLengthUpdateDisplay.test_two_atoms_enables_apply_and_shows_distance
_update_display enables Apply and shows the current distance when selection is complete._

- assert dlg.apply_button.isEnabled()
- assert 'Å' in dlg.distance_label.text()

### TestBondLengthApplyChanges.test_incomplete_selection_is_noop
_apply_changes does not raise when the atom selection is incomplete._


### TestBondLengthApplyChanges.test_invalid_input_shows_warning
_apply_changes shows a warning when the distance input is non-numeric._

- mb.warning.assert_called_once()

### TestBondLengthApplyChanges.test_negative_distance_shows_warning
_apply_changes shows a warning when a negative distance is entered._

- mb.warning.assert_called_once()

### TestBondLengthApplyChanges.test_valid_apply_pushes_undo
_apply_changes pushes an undo state when the input is valid._

- mw.edit_actions_manager.push_undo_state.assert_called()

### TestBondLengthGeometry.test_atom_only_mode_moves_atom2
_In atom-only mode atom1 stays fixed and atom2 is moved to the target distance._

- assert after_pos1 == pytest.approx(before_pos1, abs=0.0001)
- assert dist == pytest.approx(2.0, abs=0.0001)

### TestBondLengthGeometry.test_default_mode_moves_group
_In group mode atom1 stays fixed and atom2's whole side moves to the target distance._

- assert after_pos1 == pytest.approx(before_pos1, abs=0.0001)
- assert dist == pytest.approx(2.0, abs=0.0001)

### TestBondLengthGeometry.test_both_groups_mode
_In both-groups mode both halves move symmetrically to reach the target distance._

- assert dist == pytest.approx(2.0, abs=0.001)

### TestBondLengthGeometry.test_on_distance_input_changed_syncs_slider
_on_distance_input_changed updates the slider position to match the typed value._

- assert dlg.distance_slider.value() == 200

### TestAngleDialogPicking.test_sequential_picking
_Three sequential picks fill atom1, atom2, atom3 in order._

- assert dlg.atom1_idx == 2 and dlg.atom2_idx is None
- assert dlg.atom2_idx == 0 and dlg.atom3_idx is None
- assert dlg.atom3_idx == 1

### TestAngleDialogPicking.test_fourth_pick_resets
_A fourth pick after a complete selection restarts from atom1._

- assert dlg.atom1_idx == 5
- assert dlg.atom2_idx is None
- assert dlg.atom3_idx is None

### TestAngleDialogPicking.test_clear_selection
_clear_selection resets all indices and the snapshot positions to None._

- assert dlg.atom1_idx is None
- assert dlg._snapshot_positions is None

### TestAngleDialogPicking.test_preselected_atoms
_AngleDialog pre-populates all three atom indices from preselected_atoms._

- assert dlg.atom1_idx == 2
- assert dlg.atom2_idx == 0
- assert dlg.atom3_idx == 1

### TestAngleDialogIsComplete.test_incomplete_with_two_atoms
__is_selection_complete returns False when only two atoms are picked._

- assert not dlg._is_selection_complete()

### TestAngleDialogIsComplete.test_complete_with_three_atoms
__is_selection_complete returns True when all three atom indices are set._

- assert dlg._is_selection_complete()

### TestAngleDialogUpdateDisplay.test_no_atoms_disables_apply
_update_display disables the Apply button when no atoms are selected._

- assert not dlg.apply_button.isEnabled()

### TestAngleDialogUpdateDisplay.test_three_atoms_enables_apply
_update_display enables Apply and shows the angle in degrees when complete._

- assert dlg.apply_button.isEnabled()
- assert '°' in dlg.angle_label.text()

### TestAngleDialogApplyChanges.test_invalid_input_shows_warning
_apply_changes shows a warning when the angle input is non-numeric._

- mb.warning.assert_called_once()

### TestAngleDialogApplyChanges.test_angle_wraps_over_180
_270° input must wrap to -90°._

- assert float(dlg.angle_input.text()) == pytest.approx(-90.0, abs=0.01)

### TestAngleDialogApplyChanges.test_apply_changes_pushes_undo
_apply_changes pushes an undo state when the input is valid._

- mw.edit_actions_manager.push_undo_state.assert_called()

### TestAngleDialogApplyChanges.test_on_angle_input_changed_syncs_slider
_on_angle_input_changed updates the slider to match the typed angle._

- assert dlg.angle_slider.value() == 90

### TestAngleDialogGeometry.test_rotate_atom_only_mode
_In atom-only mode apply_geometry_update rotates atom3 to the target angle._

- assert new_angle == pytest.approx(120.0, abs=0.5)

### TestDihedralDialogPicking.test_sequential_picking
_Four sequential picks fill all four atom indices in order._

- assert getattr(dlg, attr) == expected

### TestDihedralDialogPicking.test_fifth_pick_resets
_A fifth pick after four atoms restarts selection with the new atom as atom1._

- assert dlg.atom1_idx == 3
- assert dlg.atom2_idx is None
- assert dlg.atom4_idx is None

### TestDihedralDialogPicking.test_clear_selection
_clear_selection resets all four atom indices and snapshot positions to None._

- assert all((getattr(dlg, a) is None for a in ['atom1_idx', 'atom2_idx', 'atom3_idx', 'atom4_idx']))
- assert dlg._snapshot_positions is None

### TestDihedralDialogPicking.test_preselected_atoms
_DihedralDialog pre-populates all four atom indices from preselected_atoms._

- assert dlg.atom1_idx == 2
- assert dlg.atom4_idx == 5

### TestDihedralDialogIsComplete.test_incomplete_with_three_atoms
__is_selection_complete returns False when only three atoms are picked._

- assert not dlg._is_selection_complete()

### TestDihedralDialogIsComplete.test_complete_with_four_atoms
__is_selection_complete returns True when all four atom indices are set._

- assert dlg._is_selection_complete()

### TestDihedralDialogCalculate.test_calculate_dihedral_incomplete_returns_zero
_calculate_dihedral returns 0.0 when the atom selection is incomplete._

- assert dlg.calculate_dihedral() == pytest.approx(0.0)

### TestDihedralDialogCalculate.test_calculate_dihedral_complete_returns_value
_calculate_dihedral returns a value in [-180, 180] for a complete selection._

- assert -180.0 <= val <= 180.0

### TestDihedralDialogUpdateDisplay.test_no_atoms_disables_apply
_update_display disables the Apply button when no atoms are selected._

- assert not dlg.apply_button.isEnabled()

### TestDihedralDialogUpdateDisplay.test_four_atoms_enables_apply
_update_display enables Apply and shows the dihedral in degrees when complete._

- assert dlg.apply_button.isEnabled()
- assert '°' in dlg.dihedral_label.text()

### TestDihedralDialogApplyChanges.test_invalid_input_shows_warning
_apply_changes shows a warning when the dihedral input is non-numeric._

- mb.warning.assert_called_once()

### TestDihedralDialogApplyChanges.test_dihedral_wraps_over_180
_270° input must wrap to -90°._

- assert float(dlg.dihedral_input.text()) == pytest.approx(-90.0, abs=0.01)

### TestDihedralDialogApplyChanges.test_apply_changes_pushes_undo
_apply_changes pushes an undo state when the dihedral input is valid._

- mw.edit_actions_manager.push_undo_state.assert_called()

### TestDihedralDialogApplyChanges.test_on_dihedral_input_changed_syncs_slider
_on_dihedral_input_changed updates the slider to match the typed dihedral._

- assert dlg.dihedral_slider.value() == 60

### TestDihedralDialogGeometry.test_rotate_atom_only_sets_dihedral
_In atom-only mode apply_geometry_update rotates atom4 to the target dihedral._

- assert result == pytest.approx(60.0, abs=1.0)

### TestDihedralDialogGeometry.test_default_group_mode_sets_dihedral
_In group mode apply_geometry_update rotates the whole side to the target dihedral._

- assert result == pytest.approx(60.0, abs=1.0)

## tests/unit/test_hydrogen.py

### test_add_hydrogen_atoms_app_logic
_Verify add_hydrogen_atoms creates items in the scene based on RDKit logic._

- assert len(h_calls) == 4
- assert actions.scene.create_bond.call_count == 4

### test_remove_hydrogen_atoms_app_logic
_Verify remove_hydrogen_atoms finds and deletes H items using app logic._

- actions.scene.delete_items.assert_called()
- assert h_item in deleted_set
- assert actions.scene.atom_items[c_id] not in deleted_set

## tests/unit/test_io.py

### test_project_save_load_logic
_Verify project save and load logic for .pmeprj files._

- assert os.path.exists(project_file)
- assert len(io_handler.data.atoms) == 1
- assert next(iter(io_handler.data.atoms.values()))['charge'] == 1

### test_xyz_export_logic
_Verify XYZ export logic._

- assert os.path.exists(xyz_file)
- assert xyz_mol.GetNumAtoms() == mol.GetNumAtoms()

## tests/unit/test_io_manager.py

### TestPromptForCharge.test_accept_returns_charge
_Accepted dialog returns the entered charge, ok=True, skip=False._

- assert charge == expected_charge
- assert ok is True
- assert skip is False

### TestPromptForCharge.test_cancel_returns_none_false_false
_Cancelled dialog returns charge=None, ok=False, skip=False._

- assert charge is None
- assert ok is False
- assert skip is False

### TestPromptForCharge.test_skip_chemistry_returns_zero_true_true
_Skip button returns charge=0, ok=True, skip=True._

- assert charge == 0
- assert ok is True
- assert skip is True

### TestLoadXYZFor3DViewing.test_explicit_path_skips_dialog
_Passing file_path directly bypasses the file dialog._

- io.load_xyz_file.assert_called_once_with(str(xyz))
- host.edit_actions_manager.clear_all.assert_called_once_with(skip_check=True)
- assert host.view_3d_manager.current_mol is mol

### TestLoadXYZFor3DViewing.test_dialog_provides_path
_Without explicit path, QFileDialog is shown and its result used._

- io.load_xyz_file.assert_called_once_with(str(xyz))
- assert host.view_3d_manager.current_mol is mol

### TestLoadXYZFor3DViewing.test_dialog_cancelled_is_noop
_Empty string from dialog → early return, no state change._

- io.load_xyz_file.assert_not_called()
- assert host.view_3d_manager.current_mol is None

### TestLoadXYZFor3DViewing.test_load_failure_shows_error
_load_xyz_file returning None shows error on status bar._

- host.statusBar_mock.showMessage.assert_called()
- assert host.view_3d_manager.current_mol is None

### TestLoadXYZFor3DViewing.test_ui_modes_enabled_on_success
_3D viewer UI mode methods are called after successful load._

- host.ui_manager.enter_3d_viewer_mode.assert_called_once()
- host.ui_manager.enable_3d_features.assert_called_once_with(True)

### TestLoadXYZFor3DViewing.test_is_xyz_derived_with_skip_prop
_is_xyz_derived=True when _xyz_skip_checks property is set on mol._

- assert host.is_xyz_derived is True

### TestLoadXYZFor3DViewing.test_is_xyz_derived_from_zero_bonds
_is_xyz_derived=True when molecule has 0 bonds (no skip flag)._

- assert host.is_xyz_derived is True

### TestLoadXYZFor3DViewing.test_current_file_path_updated
_init_manager.current_file_path is updated to the loaded file._

- assert host.init_manager.current_file_path == str(xyz)

### TestSave3DAsMol.test_no_mol_shows_error_no_dialog
_current_mol is None → error message, dialog never opened._

- mock_dlg.assert_not_called()
- host.statusBar_mock.showMessage.assert_called_with('Error: No 3D structure to save.')

### TestSave3DAsMol.test_dialog_cancelled_writes_nothing
_Empty path from dialog → no file written._

- assert not out.exists()

### TestSave3DAsMol.test_mol_file_written
_MOL file is created and non-empty._

- assert out.exists()
- assert out.stat().st_size > 0

### TestSave3DAsMol.test_header_line_replaced
_Second line is replaced with 'MoleditPy Ver. ... 3D'._

- assert len(lines) > 1
- assert 'MoleditPy Ver.' in lines[1]
- assert '3D' in lines[1]

### TestSave3DAsMol.test_success_status_message
_Status bar shows success message after saving._

- host.statusBar_mock.showMessage.assert_called()
- assert 'saved' in msg.lower() or '3D' in msg

### TestSave3DAsMol.test_exception_shows_error_message
_RuntimeError during MolToMolBlock shows error on status bar._

- host.statusBar_mock.showMessage.assert_called()
- assert 'error' in msg.lower()

## tests/unit/test_items_visual.py

### test_atom_item_visual_states
_Test AtomItem paint and boundingRect with different states._

- assert painter.drawEllipse.called
- assert painter.drawText.called
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

### test_atom_item_paint_resilience_to_deleted_bond
_Test AtomItem paint doesn't crash when a C++ bond object is deleted._

- assert len(painter.method_calls) > 0
- assert atom.bonds == [deleted_bond]

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
- assert rect.height() >= 1.0

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

- assert line.length() == 0.0

## tests/unit/test_main_qt_handler.py

### test_downgraded_pattern_routes_to_debug
_Messages matching a downgraded pattern are logged at DEBUG, not WARNING._

- mock_debug.assert_called_once_with('Qt: %s', msg)
- mock_log.assert_not_called()

### test_known_mode_maps_to_correct_level
_Each Qt message type is routed to the corresponding Python logging level._

- mock_log.assert_called_once_with(expected_level, 'Qt: %s', msg)
- mock_debug.assert_not_called()

### test_unknown_mode_falls_back_to_warning
_Unknown Qt message type defaults to logging.WARNING._

- mock_log.assert_called_once_with(logging.WARNING, 'Qt: %s', msg)

### test_read_startup_log_settings_both_true
_Returns (True, True) when both keys are set in settings.json._

- assert result == (True, True)

### test_read_startup_log_settings_defaults_false
_Returns (False, False) when keys are absent from settings.json._

- assert result == (False, False)

### test_read_startup_log_settings_missing_file
_Returns (False, False) when settings.json does not exist._

- assert result == (False, False)

### test_read_startup_log_settings_invalid_json
_Returns (False, False) on malformed JSON._

- assert result == (False, False)

## tests/unit/test_main_window_export_extended.py

### test_export_stl_success
_Test export_stl success path._

- mock_combined_mesh.save.assert_called_once_with('/path/to/export.stl', binary=True)
- window.statusBar().showMessage.assert_any_call('STL exported to /path/to/export.stl')

### test_export_stl_cancel
_Test export_stl cancellation._

- assert not any(('STL exported' in str(args) for args, _ in call_args_list))

### test_export_stl_without_current_mol_reaches_dialog
_Test export_stl does not require current_mol before choosing a file._

- assert mock_file_dialog.getSaveFileName.called
- assert not any((old_error in str(args) for args, _ in window.statusBar().showMessage.call_args_list))

### test_export_obj_mtl_success
_Test export_obj_mtl success path._

- mock_create.assert_called_once_with(meshes, '/path/to/export.obj', '/path/to/export.mtl')
- window.statusBar().showMessage.assert_any_call('OBJ+MTL files with individual colors exported to /path/to/export.obj and /path/to/export.mtl')

### test_export_color_stl_success
_Test export_color_stl success path._

- mock_combined_mesh.save.assert_called_once_with('/path/to/color_export.stl', binary=True)
- window.statusBar().showMessage.assert_any_call('STL exported to /path/to/color_export.stl')

### test_create_multi_material_obj_logic
_Test the file writing logic of create_multi_material_obj._

- mock_file.assert_any_call(obj_path, 'w')
- mock_file.assert_any_call(mtl_path, 'w')
- assert 'mtllib test.mtl' in content
- assert 'usemtl material_0_test_mesh' in content
- assert 'v 0.000000 0.000000 0.000000' in content
- assert 'f 1 2 3' in content

### test_export_from_3d_view_logic
_Test export_from_3d_view to ensure it iterates actors and extracts meshes._

- assert result is not None
- assert result.n_points == 10

### test_export_2d_png_success
_Test export_2d_png success path._

- mock_image.save.assert_called_with('/path/to/image.png', 'PNG')
- window.statusBar().showMessage.assert_any_call('2D view exported to /path/to/image.png')

### test_export_2d_svg_success
_Test export_2d_svg success path._

- mock_svg.setFileName.assert_called_with('/path/to/image.svg')
- window.statusBar().showMessage.assert_any_call('2D view exported to /path/to/image.svg')

### test_export_from_3d_view_no_color_logic
_Test the logic of extracting mesh without colors._

- assert result is not None
- assert result.n_points == 10

### test_export_from_3d_view_with_colors_logic
_Test logic of extracting mesh with colors and splitting._

- assert isinstance(results, list)
- assert len(results) >= 1

### test_create_multi_material_obj_complex_cells
_Test create_multi_material_obj with Triangle Strips and Quads._

- assert 'f 5 6 7 8' in content
- assert 'f 1 2 3' in content
- assert 'f 3 2 4' in content

## tests/unit/test_main_window_init_coverage.py

### test_imports_mainwindow
_Ensure MainWindow and its init submodule can be imported without crashing._

- assert hasattr(MainInitManager, 'init_ui')

### test_mainwindow_init_with_mocks
_Verify MainWindow instantiates MainInitManager during initialization._

- MockInitManager.assert_called_once()
- assert mw.init_manager == MockInitManager.return_value
- assert args[0] == mw

## tests/unit/test_main_window_proxies.py

### test_mainwindow_all_managers_assigned
_All manager attributes are assigned after MainWindow.__init__._

- assert hasattr(mw, 'export_manager')
- assert hasattr(mw, 'view_3d_manager')
- assert hasattr(mw, 'edit_3d_manager')
- assert hasattr(mw, 'edit_actions_manager')
- assert hasattr(mw, 'compute_manager')
- assert hasattr(mw, 'dialog_manager')
- assert hasattr(mw, 'io_manager')
- assert hasattr(mw, 'state_manager')
- assert hasattr(mw, 'string_importer_manager')
- assert hasattr(mw, 'ui_manager')
- assert hasattr(mw, 'init_manager')

### test_mainwindow_is_restoring_state_default
_is_restoring_state starts False after __init__._

- assert mw.is_restoring_state is False

### test_mainwindow_start_calculation_signal_exists
_MainWindow exposes a start_calculation signal._

- assert hasattr(MainWindow, 'start_calculation')

### test_current_mol_getter_delegates_to_view_3d_manager
_current_mol getter returns view_3d_manager.current_mol._

- assert mw.current_mol is mock_mol

### test_current_mol_setter_delegates_to_view_3d_manager
_current_mol setter stores the value on view_3d_manager._

- assert mw.view_3d_manager.current_mol is mock_mol

### test_plotter_property_delegates_to_view_3d_manager
_plotter property proxies to view_3d_manager.plotter._

- assert mw.plotter is mock_plotter

### test_data_property_delegates_to_state_manager
_data property proxies to state_manager.data._

- assert mw.data is mock_data

### test_scene_property_delegates_to_init_manager
_scene property proxies to init_manager.scene._

- assert mw.scene is mock_scene

### test_draw_molecule_3d_proxy
_draw_molecule_3d sets current_mol and delegates to view_3d_manager._

- assert mw.view_3d_manager.current_mol is mock_mol
- mw.view_3d_manager.draw_molecule_3d.assert_called_once_with(mock_mol)

### test_draw_molecule_3d_none_mol
_draw_molecule_3d with None sets current_mol to None and delegates._

- assert mw.view_3d_manager.current_mol is None
- mw.view_3d_manager.draw_molecule_3d.assert_called_once_with(None)

## tests/unit/test_measurement_calcs.py

### test_angle_dialog_logic_matches_rdkit
_Verify AngleDialog's native calculation against RDKit._

- assert app_calculated_angle == pytest.approx(rdkit_ref_angle, abs=0.01)

### test_adjust_bond_angle_internal_calc_matches_rdkit
_mol_geometry.adjust_bond_angle internally calculates the current angle before rotation._

- assert app_calculated_angle == pytest.approx(rdkit_ref_angle, abs=0.01)

### test_main_window_edit_3d_angle_logic_matches_rdkit
_Verify MainWindowEdit3d._calculate_angle matches RDKit._

- assert app_calculated_angle == pytest.approx(rdkit_ref_angle, abs=0.01)

### test_dihedral_logic_matches_rdkit
_Verify calculate_dihedral native calculation against RDKit._

- assert app_calculated_dihedral == pytest.approx(rdkit_ref_dihedral, abs=0.05)
- assert float(dialog_calculated_text) == pytest.approx(rdkit_ref_dihedral, abs=0.05)

### test_bond_length_dialog_logic_matches_rdkit
_Verify BondLengthDialog's distance calculation against RDKit._

- assert manual_distance == pytest.approx(rdkit_ref_dist, abs=0.0001)
- assert float(dialog_calculated_text) == pytest.approx(rdkit_ref_dist, abs=0.01)

### test_main_window_edit_3d_distance_logic_matches_rdkit
_Verify MainWindowEdit3d._calculate_distance matches RDKit._

- assert app_calculated_dist == pytest.approx(rdkit_ref_dist, abs=0.0001)

### test_custom_interactor_style_distance_logic_matches_rdkit
_Verify the inline distance calc logic in CustomInteractorStyle matches RDKit._

- assert app_calculated_dist == pytest.approx(rdkit_ref_dist, abs=0.0001)

### test_alignment_dialog_logic
_Verify AlignmentDialog properly aligns given atoms to the target axis._

- assert np.allclose(pos0, [0, 0, 0], atol=1e-05)
- assert pos1[0] > 0
- assert pos1[1] == pytest.approx(0.0, abs=1e-05)
- assert pos1[2] == pytest.approx(0.0, abs=1e-05)

### test_align_plane_dialog_logic
_Verify AlignPlaneDialog properly aligns selected atoms to the target plane._

- assert pos0[2] == pytest.approx(pos1[2], abs=1e-05)
- assert pos1[2] == pytest.approx(pos2[2], abs=1e-05)
- assert pos0[2] == pytest.approx(0.0, abs=1e-05)

## tests/unit/test_molecular_data.py

### test_add_atom_returns_incrementing_ids
_Verify atom IDs are sequential integers starting from 0._

- assert (id0, id1, id2) == (0, 1, 2)
- assert len(data.atoms) == 3

### test_add_atom_stores_properties
_Verify all atom properties (symbol, position, charge, radical) are stored._

- assert atom['symbol'] == 'N'
- assert atom['pos'][0] == pytest.approx(5.5)
- assert atom['pos'][1] == pytest.approx(-3.2)
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

### test_remove_atom_non_existent
_remove_atom gracefully handles removing an atom that doesn't exist._

- assert len(data.atoms) == 1

### test_remove_bond_non_existent
_remove_bond gracefully handles removing a bond that doesn't exist._

- assert len(data.bonds) == 0

### test_to_mol_block_handles_sanitization_failure
_to_mol_block falls back if RDKit molecule generation fails._

- assert mol_block is not None
- assert 'MoleditPy' in mol_block
- assert 'V2000' in mol_block

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

### test_to_template_dict
_Verify template serialization dictionary format and content._

- assert tmpl['format'] == 'PME Template'
- assert tmpl['version'] == '2.0'
- assert tmpl['application_version'] == '1.2.3'
- assert tmpl['name'] == 'Test Template'
- assert 'created' in tmpl
- assert len(tmpl['atoms']) == 2
- assert tmpl['atoms'][0]['symbol'] == 'C'
- assert tmpl['atoms'][0]['x'] == 1.0
- assert tmpl['atoms'][0]['y'] == 2.0
- assert tmpl['atoms'][0]['charge'] == 1
- assert len(tmpl['bonds']) == 1
- assert tmpl['bonds'][0]['order'] == 1
- assert tmpl['bonds'][0]['stereo'] == 1

### test_to_rdkit_mol_stereo_wedge_dash
_Verify wedge/dash stereo bonds map to BEGINWEDGE/BEGINDASH in RDKit._

- assert mol is not None
- assert wedge_bond.GetBondDir() == Chem.BondDir.BEGINWEDGE
- assert dash_bond.GetBondDir() == Chem.BondDir.BEGINDASH

### test_to_rdkit_mol_ez_stereo
_Verify E/Z stereo double bond maps to STEREOZ in RDKit._

- assert mol is not None
- assert double_bond.GetBondType() == Chem.BondType.DOUBLE
- assert double_bond.GetStereo() == Chem.BondStereo.STEREOZ

## tests/unit/test_molecule_scene_behavior.py

### TestGetSetting.test_returns_value_from_settings
_get_setting returns the value stored in init_manager.settings._

- assert scene.get_setting('my_key') == 'my_value'

### TestGetSetting.test_returns_default_when_key_missing
_get_setting returns the default value when the key is absent._

- assert scene.get_setting('no_such_key', 'fallback') == 'fallback'

### TestGetSetting.test_returns_default_when_window_is_none
_get_setting returns the default value when scene.window is None._

- assert scene.get_setting('any_key', 99) == 99

### TestUpdateConnectedBonds.test_calls_update_position_on_bond
_update_connected_bonds calls update_position on each live bond._

- bond.update_position.assert_called_once()

### TestUpdateConnectedBonds.test_skips_sip_deleted_bonds
_update_connected_bonds skips bonds that have been sip-deleted._

- bond.update_position.assert_not_called()

### TestUpdateConnectedBonds.test_deduplicates_shared_bond
_Bond shared by two atoms should only be updated once._

- bond.update_position.assert_called_once()

### TestUpdateConnectedBonds.test_handles_atom_without_bonds_attribute
_update_connected_bonds does not raise when an atom has no bonds attribute._


### TestClearAllProblemFlags.test_returns_true_when_flags_were_set
_clear_all_problem_flags returns True when at least one atom had has_problem set._

- assert scene.clear_all_problem_flags() is True

### TestClearAllProblemFlags.test_returns_false_when_no_flags_set
_clear_all_problem_flags returns False when no atom has has_problem set._

- assert scene.clear_all_problem_flags() is False

### TestClearAllProblemFlags.test_flag_is_reset_to_false
_clear_all_problem_flags resets has_problem to False on each atom._

- assert item.has_problem is False

### TestPurgeDeletedItems.test_noop_on_empty_list
_purge_deleted_items is a no-op when the deleted items list is empty._


### TestPurgeDeletedItems.test_clears_the_list
_purge_deleted_items empties the _deleted_items list._

- assert scene._deleted_items == []

### TestPurgeDeletedItems.test_calls_hide_on_valid_objects
_purge_deleted_items calls hide() on each live item._

- obj.hide.assert_called_once()

### TestPurgeDeletedItems.test_skips_already_sip_deleted
_purge_deleted_items skips items that have already been sip-deleted._

- obj.hide.assert_not_called()

### TestSetHoveredItem.test_stores_item
_set_hovered_item stores the provided item reference._

- assert scene.hovered_item is item

### TestSetHoveredItem.test_accepts_none
_set_hovered_item accepts None to clear the hovered item._

- assert scene.hovered_item is None

### TestEZStereoCycling.test_none_to_z_on_first_click
_Clicking a double bond with no stereo sets it to Z (stereo=3)._

- mock_stereo.assert_called_once_with(bond_item, 3)

### TestEZStereoCycling.test_z_to_e_on_click
_Clicking a Z double bond advances it to E (stereo=4)._

- mock_stereo.assert_called_once_with(bond_item, 4)

### TestEZStereoCycling.test_e_to_none_on_click
_Clicking an E double bond cycles it back to no stereo (stereo=0)._

- mock_stereo.assert_called_once_with(bond_item, 0)

### TestBondDirectionInversion.test_stereo_bond_click_inverts_atom_order
_Clicking a stereo bond with matching stereo type inverts the atom order._

- assert (a2_id, a1_id) in mock_parser_host.state_manager.data.bonds

### TestDoubleClickSelectMode.test_selects_connected_atom
_Double-click in select mode selects the connected component via BFS._

- atom_item.setSelected.assert_called_with(True)

### TestDoubleClickSelectMode.test_bond_2_5_mode_accepts_event
_Double-click in bond_2_5 mode accepts the event and does not raise._

- ev.accept.assert_called()

### TestDoubleClickChargeRadical.test_radical_increments_on_double_click
_Double-clicking an atom in radical mode increments its radical count._

- assert item.radical == 1

### TestDoubleClickChargeRadical.test_charge_plus_increments_on_double_click
_Double-clicking an atom in charge_plus mode increments its charge._

- assert item.charge == 1

### TestDoubleClickChargeRadical.test_charge_minus_decrements_on_double_click
_Double-clicking an atom in charge_minus mode decrements its charge._

- assert item.charge == -1

## tests/unit/test_move_group_dialog.py

### TestOnAtomPicked.test_picks_entire_connected_component
_Picking any atom in ethane should BFS to all 8 atoms._

- assert len(dlg.group_atoms) == mol.GetNumAtoms()

### TestOnAtomPicked.test_picking_again_toggles_deselect
_Re-picking any atom in the already-selected group deselects everything._

- assert len(dlg.group_atoms) == 0

### TestOnAtomPicked.test_bfs_stays_within_fragment
_In a two-fragment system, BFS must not cross to the other fragment._

- assert len(frags) == 2
- assert len(dlg.group_atoms) == len(frag_a)
- assert dlg.group_atoms == frag_a
- assert dlg.group_atoms.isdisjoint(frag_b)

### TestOnAtomPicked.test_selected_atoms_records_clicked_atom
_on_atom_picked records the clicked atom index in selected_atoms._

- assert 3 in dlg.selected_atoms

### TestOnAtomPicked.test_skip_pick_when_dragging
_on_atom_picked must be a no-op while dragging._

- mock_labels.assert_not_called()
- assert len(dlg.group_atoms) == 0

### TestUpdateDisplay.test_no_group_shows_placeholder
_update_display shows 'No group' when no atoms are selected._

- assert 'No group' in dlg.selection_label.text()

### TestUpdateDisplay.test_group_shows_count_and_symbols
_update_display shows atom count and symbols when a group is selected._

- assert '8' in text or 'atoms' in text.lower()

### TestUpdateDisplay.test_more_than_5_atoms_appended_with_ellipsis
_If >5 atoms selected, display must show '...' at the end._

- assert '...' in dlg.selection_label.text()

### TestApplyTranslation.test_no_group_shows_warning
_apply_translation shows a warning when no group is selected._

- mb.warning.assert_called_once()

### TestApplyTranslation.test_invalid_input_shows_warning
_apply_translation shows a warning when a translation input is non-numeric._

- mb.warning.assert_called_once()

### TestApplyTranslation.test_translation_updates_conformer_positions
_apply_translation shifts all group atoms by the given XYZ vector._

- assert after[idx] == pytest.approx(before[idx] + [1.0, 2.0, 3.0], abs=0.0001)

### TestApplyTranslation.test_translation_pushes_undo
_apply_translation pushes an undo state after moving the group._

- mw.edit_actions_manager.push_undo_state.assert_called()

### TestApplyRotation.test_no_group_shows_warning
_apply_rotation shows a warning when no group is selected._

- mb.warning.assert_called_once()

### TestApplyRotation.test_invalid_input_shows_warning
_apply_rotation shows a warning when a rotation input is non-numeric._

- mb.warning.assert_called_once()

### TestApplyRotation.test_zero_rotation_leaves_positions_unchanged
_apply_rotation with all-zero angles does not move any atom._

- assert after == pytest.approx(before, abs=1e-05)

### TestApplyRotation.test_90deg_z_rotation_around_centroid
_90° rotation around Z maps (centroid+[1,0,0]) → (centroid+[0,1,0])._

- assert after[0] == pytest.approx([0.0, 1.0, 0.0], abs=1e-05)
- assert after[1] == pytest.approx([0.0, -1.0, 0.0], abs=1e-05)

### TestApplyRotation.test_rotation_pushes_undo
_apply_rotation pushes an undo state after rotating the group._

- mw.edit_actions_manager.push_undo_state.assert_called()

### TestResetInputs.test_reset_translation_inputs
_reset_translation_inputs sets all three translation fields to '0.0'._

- assert dlg.x_trans_input.text() == '0.0'
- assert dlg.y_trans_input.text() == '0.0'
- assert dlg.z_trans_input.text() == '0.0'

### TestResetInputs.test_reset_rotation_inputs
_reset_rotation_inputs sets all three rotation fields to '0.0'._

- assert dlg.x_rot_input.text() == '0.0'
- assert dlg.y_rot_input.text() == '0.0'
- assert dlg.z_rot_input.text() == '0.0'

### TestClearSelection.test_clear_removes_group_and_selected_atoms
_clear_selection empties group_atoms and selected_atoms._

- assert len(dlg.group_atoms) == 0
- assert len(dlg.selected_atoms) == 0

### TestClearSelection.test_clear_resets_drag_state
_clear_selection resets is_dragging_group and drag_start_pos to their default values._

- assert not dlg.is_dragging_group
- assert dlg.drag_start_pos is None

### TestClearSelection.test_clear_updates_display
_clear_selection updates the display label back to 'No group'._

- assert 'No group' in dlg.selection_label.text()

### test_show_atom_labels_camera_restore
_show_atom_labels restores the camera position after adding highlight meshes._

- mock_clear.assert_called_once()
- assert kwargs.get('reset_camera') is False
- assert plotter.camera_position == [(1, 2, 3), (4, 5, 6), (7, 8, 9)]
- plotter.render.assert_called()

### TestInit.test_preselected_atoms_triggers_on_atom_picked
_Passing preselected_atoms to __init__ pre-selects the connected group via BFS._

- assert len(dlg.group_atoms) == mol.GetNumAtoms()

### TestInit.test_no_preselected_atoms_leaves_group_empty
_Without preselected_atoms the group_atoms set is empty on init._

- assert len(dlg.group_atoms) == 0

### TestInit.test_initial_drag_state
_MoveGroupDialog initialises with dragging flags all False/None._

- assert dlg.is_dragging_group is False
- assert dlg.drag_start_pos is None
- assert dlg.potential_drag is False

### TestInit.test_window_title
_MoveGroupDialog window title is 'Move Group'._

- assert dlg.windowTitle() == 'Move Group'

### TestInitUI.test_translation_inputs_default_zero
_Translation input fields default to '0.0' on dialog creation._

- assert dlg.x_trans_input.text() == '0.0'
- assert dlg.y_trans_input.text() == '0.0'
- assert dlg.z_trans_input.text() == '0.0'

### TestInitUI.test_rotation_inputs_default_zero
_Rotation input fields default to '0.0' on dialog creation._

- assert dlg.x_rot_input.text() == '0.0'
- assert dlg.y_rot_input.text() == '0.0'
- assert dlg.z_rot_input.text() == '0.0'

### TestInitUI.test_selection_label_initial_text
_The selection label shows 'No group' before any atom is picked._

- assert 'No group' in dlg.selection_label.text()

### TestUpdateDisplayBoundary.test_exactly_5_atoms_no_ellipsis
_Exactly 5 selected atoms must NOT show '...'._

- assert '5' in text
- assert '...' not in text

### TestUpdateDisplayBoundary.test_more_than_5_atoms_has_ellipsis
_8-atom ethane (with H) must show '...' after 5th._

- assert '...' in dlg.selection_label.text()

### TestApplyRotationAxes.test_90deg_x_rotation_around_centroid
_90° rotation around X maps [0,1,0] → [0,0,1] relative to centroid._

- assert after[0] == pytest.approx([0.0, 0.0, 1.0], abs=1e-05)
- assert after[1] == pytest.approx([0.0, 0.0, -1.0], abs=1e-05)

### TestApplyRotationAxes.test_90deg_y_rotation_around_centroid
_90° rotation around Y maps [1,0,0] → [0,0,-1] relative to centroid._

- assert after[0] == pytest.approx([0.0, 0.0, -1.0], abs=1e-05)
- assert after[1] == pytest.approx([0.0, 0.0, 1.0], abs=1e-05)

### TestApplyRotationAxes.test_combined_rotation_pushes_undo
_apply_rotation with combined XYZ angles pushes an undo state._

- mw.edit_actions_manager.push_undo_state.assert_called()

### TestAtomLabels.test_show_atom_labels_calls_plotter_add_mesh
_show_atom_labels adds highlight spheres to the plotter and calls render._

- mock_plotter.add_mesh.assert_called()
- mock_plotter.render.assert_called()

### TestAtomLabels.test_show_atom_labels_no_group_does_nothing
_show_atom_labels does not call add_mesh when group_atoms is empty._

- mw.view_3d_manager.plotter.add_mesh.assert_not_called()

### TestAtomLabels.test_clear_atom_labels_removes_highlight_actor
_clear_atom_labels removes the highlight actor from the plotter and sets it to None._

- mock_plotter.remove_actor.assert_called()
- assert dlg.highlight_actor is None

### TestAtomLabels.test_clear_atom_labels_none_plotter_does_not_raise
_clear_atom_labels does not raise when the plotter is None._


### TestEventFilter.test_returns_false_when_plotter_is_none
_eventFilter returns False immediately when the plotter is None._

- assert result is False

### TestEventFilter.test_returns_false_when_mol_is_none
_eventFilter returns False immediately when the molecule is None._

- assert result is False

### TestEventFilter.test_double_click_resets_state_and_returns_false
_A double-click resets dragging flags and returns False._

- assert result is False
- assert dlg.is_dragging_group is False
- assert dlg.potential_drag is False

### TestEventFilter.test_non_interactor_obj_delegates_to_super
_Events on objects other than the plotter interactor use base behaviour._

- assert result is False

### TestMoveGroupDeselectToggle.test_on_atom_picked_deselects_connected_group
_Re-picking the same atom deselects the entire connected group._

- assert len(dlg.group_atoms) > 0
- assert 0 in dlg.selected_atoms
- assert len(dlg.group_atoms) == 0
- assert 0 not in dlg.selected_atoms

## tests/unit/test_move_selected_atoms_dialog.py

### TestOnAtomPicked.test_picks_atom_individually_no_bfs
_Picking an atom should select ONLY that atom (no BFS connected components)._

- assert dlg.selected_atoms == {0}

### TestOnAtomPicked.test_picking_again_toggles_deselect
_Re-picking the same atom deselects it._

- assert 0 in dlg.selected_atoms
- assert 0 not in dlg.selected_atoms

### TestApplyTranslation.test_no_atoms_shows_warning
_apply_translation shows a warning when no atoms are selected._

- mb.warning.assert_called_once()

### TestApplyTranslation.test_invalid_input_shows_warning
_apply_translation shows a warning when a translation input is non-numeric._

- mb.warning.assert_called_once()

### TestApplyTranslation.test_translation_updates_only_selected_atoms
_apply_translation moves only selected atoms and leaves others unchanged._

- assert after[0] == pytest.approx(before[0] + [1.0, 2.0, 3.0], abs=0.0001)
- assert after[1] == pytest.approx(before[1], abs=0.0001)

### TestApplyRotation.test_no_atoms_shows_warning
_apply_rotation shows a warning when no atoms are selected._

- mb.warning.assert_called_once()

### TestApplyRotation.test_invalid_input_shows_warning
_apply_rotation shows a warning when a rotation input is non-numeric._

- mb.warning.assert_called_once()

### TestApplyRotation.test_rotation_updates_only_selected_atoms_around_centroid
_apply_rotation rotates only selected atoms around their centroid._

- assert after[0] == pytest.approx([0.0, 1.0, 0.0], abs=1e-05)
- assert after[1] == pytest.approx([0.0, -1.0, 0.0], abs=1e-05)
- assert after[2] == pytest.approx([5.0, 5.0, 5.0], abs=1e-05)

### TestResetInputs.test_reset_translation_inputs
_reset_translation_inputs sets all three translation fields to '0.0'._

- assert dlg.x_trans_input.text() == '0.0'
- assert dlg.y_trans_input.text() == '0.0'
- assert dlg.z_trans_input.text() == '0.0'

### TestResetInputs.test_reset_rotation_inputs
_reset_rotation_inputs sets all three rotation fields to '0.0'._

- assert dlg.x_rot_input.text() == '0.0'
- assert dlg.y_rot_input.text() == '0.0'
- assert dlg.z_rot_input.text() == '0.0'

### TestClearSelection.test_clear_removes_selected_atoms
_clear_selection empties the selected_atoms set._

- assert len(dlg.selected_atoms) == 0

### test_group_atoms_property
_group_atoms property aliases selected_atoms for compatibility with the mixin._

- assert dlg.selected_atoms == {0, 2}
- assert dlg.group_atoms == {0, 2}

### test_preselected_atoms_init
_MoveSelectedAtomsDialog pre-populates selected_atoms from preselected_atoms._

- assert dlg.selected_atoms == {0, 1}

### test_show_atom_labels_none_positions
_show_atom_labels logs an error and does nothing when atom_positions_3d is None._

- mock_log.assert_called_once_with('atom_positions_3d is None in update_atom_labels')

### test_show_and_clear_atom_labels
_show_atom_labels calls add_mesh; clear_atom_labels removes the actor._

- plotter.remove_actor.assert_called()
- plotter.add_mesh.assert_called_once()
- plotter.render.assert_called()

### test_event_filter_unrelated_obj
_eventFilter returns False for events on objects other than the VTK interactor._

- assert dlg.eventFilter(QWidget(), event) is False

### test_event_filter_double_click
_A double-click event resets potential_drag and returns False._

- assert dlg.eventFilter(plotter.interactor, event) is False
- assert dlg.potential_drag is False

### test_event_filter_mouse_press_with_selection
_A left-press with atoms already selected delegates to VTK (returns False)._

- assert dlg.eventFilter(plotter.interactor, event) is False

### test_event_filter_mouse_press_selects_atom
_A left-press that hits an atom picks it and returns True._

- assert dlg.eventFilter(plotter.interactor, event) is True
- mock_pick.assert_called_once_with(0)
- assert dlg._consume_next_left_release is True

### test_event_filter_mouse_press_empty_space
_A left-press that misses all atoms returns False._

- assert dlg.eventFilter(plotter.interactor, event) is False

### test_event_filter_mouse_move_hover_cursor
_Mouse-move over a selected atom changes the cursor to OpenHand; otherwise to Arrow._

- plotter.setCursor.assert_called_with(Qt.CursorShape.OpenHandCursor)
- plotter.setCursor.assert_called_with(Qt.CursorShape.ArrowCursor)

### test_event_filter_mouse_release_consume
_A mouse-release is consumed and clears _consume_next_left_release when the flag is set._

- assert dlg.eventFilter(plotter.interactor, event) is True
- assert dlg._consume_next_left_release is False

### test_mouse_move_potential_drag_to_actual_drag
_A mouse-move beyond the threshold converts a potential drag to an actual drag._

- assert dlg.eventFilter(plotter.interactor, event) is True
- assert dlg.is_dragging_group is True
- assert dlg.potential_drag is False
- plotter.setCursor.assert_called_with(Qt.CursorShape.ClosedHandCursor)

### test_mouse_move_during_actual_drag
_Mouse-move while is_dragging_group is True sets mouse_moved_during_drag._

- assert dlg.eventFilter(plotter.interactor, event) is True
- assert dlg.mouse_moved_during_drag is True

### test_mouse_release_no_movement_toggles_atom
_A release with no movement calls on_atom_picked to toggle the atom selection._

- assert dlg.eventFilter(plotter.interactor, event) is True
- mock_pick.assert_called_once_with(0)
- assert dlg.potential_drag is False

### test_mouse_release_with_movement_resets_drag_state
_A release after dragging resets drag state without calling on_atom_picked._

- assert dlg.eventFilter(plotter.interactor, event) is True
- mock_pick.assert_not_called()
- assert dlg.is_dragging_group is False

### test_handle_mouse_press_exceptions
_eventFilter returns False and does not raise when GetEventPosition throws._

- assert dlg.eventFilter(plotter.interactor, event) is False

### test_on_atom_picked_during_drag_ignored
_on_atom_picked is a no-op while is_dragging_group is True._

- assert 0 not in dlg.selected_atoms

### test_update_display_many_atoms
_update_display shows '...' and total count when more than 5 atoms are selected._

- assert '...' in text
- assert 'Selected: 7 atoms' in text

### TestClickToDeselect.test_eventfilter_delegates_to_vtk_when_atoms_selected
_eventFilter returns False to let VTK handle the press when atoms are selected._

- assert result is False

### TestClickToDeselect.test_eventfilter_handles_press_when_no_atoms_selected
_A left-press with no prior selection picks the clicked atom and returns True._

- assert result is True
- assert 0 in dlg.selected_atoms

### TestClickToDeselect.test_eventfilter_release_click_only_deselects_atom
_A release after a pure click (no drag) calls on_atom_picked to toggle the atom._

- assert result is True
- mock_on_pick.assert_called_once_with(0)

### test_show_atom_labels_camera_restore
_show_atom_labels restores the camera position after adding highlight meshes._

- mock_clear.assert_called_once()
- assert kwargs.get('reset_camera') is False
- assert plotter.camera_position == [(1, 2, 3), (4, 5, 6), (7, 8, 9)]
- plotter.render.assert_called()

## tests/unit/test_parser_robustness.py

### test_set_mol_prop_safe_robustness
_Test that _set_mol_prop handles RDKit exceptions by narrowing._

- assert not mol.HasProp('test_prop')

### test_get_set_mol_prop_instance
_Test instance versions of prop helpers._

- assert win._get_mol_prop(mol, 'int_prop') == 123
- assert pytest.approx(win._get_mol_prop(mol, 'float_prop')) == 1.23
- assert win._get_mol_prop(mol, 'str_prop') == 'hello'

### test_save_as_xyz_no_mol
_Test save_as_xyz handles missing molecule gracefully._

- win.statusBar().showMessage.assert_called_with('Error: Please generate a 3D structure first.')

### test_save_as_xyz_logic
_Test the core logic of save_as_xyz (mocking the file dialog)._

- assert os.path.exists(save_path)
- assert '1' in content
- assert 'Generated by MoleditPy' in content

### test_load_mol_counts_fix
_Test that the MOL counts line fixer is robust._

- assert 'V2000' in fixed
- assert len(fixed) == 39
- assert win.fix_mol_counts_line(line_v2) == line_v2

### test_load_xyz_robustness
_Test that XYZ parser handles malformed files via statusBar._

- win.statusBar().showMessage.assert_called()
- assert 'Error parsing XYZ file' in args[0]

## tests/unit/test_parsers.py

### test_fix_mol_block
_Verify the logic that fixes malformed MOL block counts lines._

- assert 'V2000' in lines[3]
- assert len(lines[3]) >= 39

### test_load_mol_file_logic
_Verify loading of a standard MOL file._

- assert len(parser.data.atoms) == 2
- assert any((a['symbol'] == 'C' for a in parser.data.atoms.values()))

### test_xyz_parsing_logic
_Verify basic XYZ parsing logic._

- assert mol is not None
- assert mol.GetNumAtoms() == 3
- assert any((a.GetSymbol() == 'O' for a in mol.GetAtoms()))

### test_load_xyz_file_with_estimation
_Verify XYZ loading with automatic bond estimation._

- assert mol is not None
- assert mol.GetNumAtoms() == 2
- assert mol.GetNumBonds() == 1

### test_save_as_mol_logic
_Verify saving a molecule as a MOL file._

- assert os.path.exists(save_path)
- assert 'C  ' in content

### test_load_mol_file_fallback_to_sd_supplier
_Verify fallback to ForwardSDMolSupplier when standard MolBlock reading fails._

- assert len(parser.data.atoms) == 1

### test_load_xyz_always_ask_charge
_Verify that charge dialog is shown when loading XYZ._

- assert mol is not None
- assert mol.GetIntProp('_xyz_charge') == 1

### test_load_xyz_charge_loop_cancel
_Verify handling of user cancellation during the charge input loop._

- assert result is None

### test_load_xyz_unrecognized_symbol
_Test load_xyz_file raises ValueError for unrecognized element symbols._

- assert mol is None
- assert any(('Unrecognized element symbol' in m for m in msgs))

### test_load_mol_file_with_v2000_fix
_Verify that malformed V2000 headers are fixed automatically during MOL load._

- assert len(parser.data.atoms) == 1

### test_load_xyz_recovery_loop_retries
_Verify that the XYZ load recovery loop handles retries correctly._

- assert mol is not None
- assert mol.GetIntProp('_xyz_charge') == 1

### test_load_mol_file_not_found
_Verify error handling when a MOL file is not found._

- parser.statusBar().showMessage.assert_any_call('File not found: missing_parser_xyz_final.mol')

### test_load_mol_file_invalid_format
_Verify error handling for invalid MOL file content._

- assert any(('Invalid MOL' in m or 'Failed to read' in m or 'Error loading' in m for m in msgs))

### test_save_as_xyz_charge_mult
_Verify that charge and multiplicity are written to XYZ file headers._

- assert 'chrg = 1' in f.read()

### test_load_mol_file_malformed_counts
_Verify fixing of malformed counts lines in MOL files._

- assert len(parser.data.atoms) == 1

### test_load_xyz_complex_recovery_branches
_Verify complex error recovery branches during XYZ loading._

- assert mol is not None
- assert mol.HasProp('_xyz_skip_checks') or getattr(mol, '_xyz_skip_checks', False)

### test_save_as_mol_no_current_path
_Verify saving as MOL when no current file path is set._

- assert os.path.exists(save_path)

### test_load_xyz_skip_chemistry_via_button
_Verify skip chemistry check flag is set when user chooses to skip._

- assert mol is not None
- assert mol.HasProp('_xyz_skip_checks') or getattr(mol, '_xyz_skip_checks', False)

## tests/unit/test_periodic_table_dialog.py

### test_periodic_table_dialog_no_parent
_Test dialog creation without parent settings._

- assert dialog.windowTitle() == 'Select an Element'
- assert dialog.layout() is not None

### test_periodic_table_dialog_with_parent_overrides
_Test dialog creation with parent settings cpk_colors overrides._

- assert dialog.windowTitle() == 'Select an Element'

### test_periodic_table_dialog_with_parent_error_handling
_Test dialog creation when parent attribute access raises exception._

- assert dialog.windowTitle() == 'Select an Element'

### test_periodic_table_dialog_element_clicked
_Test that clicking a button emits element_selected and accepts the dialog._

- assert blocker.args == ['H']
- assert h_button is not None

## tests/unit/test_planarize_dialog.py

### TestApplyPlanarizeGuard.test_fewer_than_three_atoms_shows_warning
_Warning shown when fewer than 3 atoms are selected for planarize._

- mb.warning.assert_called_once()

### TestApplyPlanarizeGuard.test_zero_atoms_shows_warning
_Warning shown when no atoms are selected for planarize._

- mb.warning.assert_called_once()

### TestApplyPlanarizeGeometry.test_planarize_reduces_z_spread
_Planarize reduces the Z-coordinate variance of selected atoms._

- assert after_z_var <= before_z_var + 1e-06

### TestApplyPlanarizeGeometry.test_planarize_pushes_undo
_Planarize pushes an undo state after geometry update._

- mw.edit_actions_manager.push_undo_state.assert_called()

### TestApplyPlanarizeGeometry.test_apply_calls_draw_molecule_3d
_Planarize triggers a 3D redraw after modifying geometry._

- mw.view_3d_manager.draw_molecule_3d.assert_called()

## tests/unit/test_plugin_context.py

### TestMarkProjectModified.test_sets_has_unsaved_changes
_mark_project_modified sets has_unsaved_changes to True on the state manager._

- self.assertTrue(mw.state_manager.has_unsaved_changes)

### TestMarkProjectModified.test_calls_update_window_title
_mark_project_modified calls update_window_title on the state manager._

- mw.state_manager.update_window_title.assert_called_once()

### TestMarkProjectModified.test_no_crash_when_state_manager_missing
_mark_project_modified does not raise when the main window has no state_manager._


### TestMarkProjectModified.test_no_crash_when_update_window_title_missing
_mark_project_modified does not raise when update_window_title is absent._


### TestMarkProjectModified.test_no_crash_when_state_manager_raises
_mark_project_modified does not propagate exceptions from state_manager._


### TestMarkProjectModified.test_method_exists_on_plugincontext
_PluginContext exposes a callable mark_project_modified method._

- self.assertTrue(hasattr(PluginContext, 'mark_project_modified'), 'PluginContext must expose mark_project_modified()')
- self.assertTrue(callable(getattr(PluginContext, 'mark_project_modified')))

### TestRefreshUi.test_calls_all_required_managers
_refresh_ui calls update_realtime_info, update_undo_redo_actions, and update_window_title._

- mw.state_manager.update_realtime_info.assert_called_once()
- mw.edit_actions_manager.update_undo_redo_actions.assert_called_once()
- mw.state_manager.update_window_title.assert_called_once()

### TestRefreshUi.test_no_crash_when_managers_missing
_refresh_ui does not raise when the main window has no manager attributes._


### TestFit3dView.test_calls_fit_to_view
_fit_3d_view delegates to view_3d_manager.fit_to_view._

- mw.view_3d_manager.fit_to_view.assert_called_once()

### TestFit3dView.test_no_crash_when_fit_to_view_missing
_fit_3d_view does not raise when fit_to_view is absent on the view manager._


### test_clear_canvas_delegates
_clear_canvas delegates to edit_actions_manager.clear_2d_editor with push_to_undo._

- mw.edit_actions_manager.clear_2d_editor.assert_called_once_with(push_to_undo=push_to_undo)

### test_clear_canvas_no_crash_when_manager_missing
_clear_canvas does not raise when edit_actions_manager is absent._


### test_set_3d_features_enabled_delegates
_set_3d_features_enabled delegates to ui_manager.enable_3d_features._

- mw.ui_manager.enable_3d_features.assert_called_once_with(enabled)

### test_set_3d_features_enabled_no_crash_when_ui_manager_missing
_set_3d_features_enabled does not raise when ui_manager is absent._


### test_set_analysis_enabled_delegates
_set_analysis_enabled delegates to init_manager.analysis_action.setEnabled._

- mw.init_manager.analysis_action.setEnabled.assert_called_once_with(enabled)

### test_set_analysis_enabled_no_crash_when_action_missing
_set_analysis_enabled does not raise when analysis_action is absent._


### TestCheckChemistryProblems.test_calls_fallback
_check_chemistry_problems delegates to compute_manager.check_chemistry_problems_fallback._

- mw.compute_manager.check_chemistry_problems_fallback.assert_called_once()

### TestCheckChemistryProblems.test_no_crash_when_compute_manager_missing
_check_chemistry_problems does not raise when compute_manager is absent._


### TestRefresh2dScene.test_calls_update_all_items
_refresh_2d_scene delegates to init_manager.scene.update_all_items._

- mw.init_manager.scene.update_all_items.assert_called_once()

### TestRefresh2dScene.test_no_crash_when_init_manager_missing
_refresh_2d_scene does not raise when init_manager is absent._


### TestRefresh2dScene.test_no_crash_when_scene_missing
_refresh_2d_scene does not raise when scene is absent on init_manager._


### TestRefresh2dScene.test_no_crash_when_update_all_items_missing
_refresh_2d_scene does not raise when update_all_items is absent on scene._


### test_no_crash_when_mw_is_none
_All PluginContext API methods are no-ops when the main window is None._


## tests/unit/test_plugin_interface.py

### TestPluginInterface.test_plugin_context_init
_Test PluginContext initialization._

- assert ctx._manager == mock_manager
- assert ctx._plugin_name == 'TestPlugin'

### TestPluginInterface.test_delegates_to_manager
_Each PluginContext API call delegates to the corresponding manager method._


### TestPluginInterface.test_get_3d_controller
_Test get_3d_controller returns a controller linked to main window._

- assert isinstance(controller, Plugin3DController)
- assert controller._mw == mock_main_window

### TestPluginInterface.test_get_main_window
_Test get_main_window delegation._

- assert ctx.get_main_window() == mock_main_window

### TestPluginInterface.test_current_molecule_property
_Test current_molecule getter and setter._

- assert ctx.current_molecule == 'mock_molecule'
- assert mock_main_window.view_3d_manager.current_mol == 'new_molecule'
- mock_main_window.view_3d_manager.draw_molecule_3d.assert_called_once_with('new_molecule')

### TestPluginInterface.test_current_molecule_no_window
_Test current_molecule when main window is None._

- assert ctx.current_molecule is None

### TestPluginInterface.test_register_3d_context_menu
_Test deprecated register_3d_context_menu emits DeprecationWarning._


### TestPluginInterface.test_3d_controller_set_atom_color
_Test Plugin3DController.set_atom_color._

- mock_main_window.view_3d_manager.update_atom_color_override.assert_called_once_with(1, '#FF0000')
- mock_main_window.view_3d_manager.plotter.render.assert_called_once()

### TestPluginInterface.test_3d_controller_set_bond_color
_Test Plugin3DController.set_bond_color._

- mock_main_window.view_3d_manager.update_bond_color_override.assert_called_once_with(2, '#00FF00')
- mock_main_window.view_3d_manager.plotter.render.assert_called_once()

### TestPluginInterface.test_add_plugin_menu
_add_plugin_menu prepends 'Plugin/' to the path._

- mock_manager.register_menu_action.assert_called_once_with('TestPlugin', 'Plugin/Utility/My Tool...', callback, None, None, None)

### TestPluginInterface.test_add_plugin_menu_strips_leading_slash
_add_plugin_menu strips a leading slash from the path._

- mock_manager.register_menu_action.assert_called_once_with('TestPlugin', 'Plugin/Analysis/Viewer', callback, None, None, None)

### TestPluginInterface.test_add_plugin_menu_with_text_and_shortcut
_add_plugin_menu passes optional text/icon/shortcut through._

- mock_manager.register_menu_action.assert_called_once_with('TestPlugin', 'Plugin/File/Export...', callback, 'Export', None, 'Ctrl+E')

### TestPluginInterface.test_register_menu_action_new_style
_register_menu_action (new style: path, callback) delegates correctly._

- mock_manager.register_menu_action.assert_called_once_with('TestPlugin', 'File/Open', callback, None, None, None)

### TestPluginInterface.test_register_menu_action_old_style
_register_menu_action (old style: path, text, callback) delegates correctly._

- mock_manager.register_menu_action.assert_called_once_with('TestPlugin', 'File/Import', callback, 'Import PubChem...', None, None)

### TestPluginInterface.test_get_setting_returns_default_when_missing
_get_setting returns the default if the key is absent._

- assert ctx.get_setting('theme', 'light') == 'light'

### TestPluginInterface.test_get_setting_returns_stored_value
_get_setting returns the stored value when the namespaced key exists._

- assert ctx.get_setting('theme', 'light') == 'dark'

### TestPluginInterface.test_get_setting_namespacing
_get_setting is namespaced — same key for different plugins is independent._

- assert ctx_a.get_setting('color') == 'red'
- assert ctx_b.get_setting('color') == 'blue'

### TestPluginInterface.test_get_setting_no_main_window
_get_setting returns default when main window is None._

- assert ctx.get_setting('key', 'fallback') == 'fallback'

### TestPluginInterface.test_set_setting_writes_namespaced_key
_set_setting writes to init_manager.settings with correct namespace._

- assert mw.init_manager.settings['plugin.MyPlugin.theme'] == 'dark'

### TestPluginInterface.test_set_setting_marks_dirty
_set_setting sets settings_dirty = True._

- assert mw.init_manager.settings_dirty is True

### TestPluginInterface.test_set_setting_overwrites_existing
_set_setting overwrites an existing value._

- assert mw.init_manager.settings['plugin.MyPlugin.x'] == 99

### TestPluginInterface.test_set_setting_no_main_window
_set_setting is a no-op when main window is None._


### TestPluginInterface.test_get_after_set_roundtrip
_Value written with set_setting can be read back with get_setting._

- assert ctx.get_setting('count', 0) == 7

### TestPluginInterface.test_plotter
_ctx.plotter proxies to view_3d_manager.plotter, or None if no window._

- assert ctx.plotter == 'mock_plotter'
- assert ctx.plotter is None

### TestPluginInterface.test_scene
_ctx.scene proxies to init_manager.scene, or None if no window._

- assert ctx.scene == 'mock_scene'
- assert ctx.scene is None

### TestPluginInterface.test_draw_molecule_3d
_draw_molecule_3d delegates to view_3d_manager and is safe with no window._

- mock_main_window.view_3d_manager.draw_molecule_3d.assert_called_once_with('mol')

### TestPluginInterface.test_refresh_3d_view
_refresh_3d_view redraws the molecule, or just renders if no mol._

- mock_main_window.view_3d_manager.draw_molecule_3d.assert_called_once_with('mol')
- mock_main_window.view_3d_manager.plotter.render.assert_called_once()

### TestPluginInterface.test_reset_3d_camera
_reset_3d_camera calls plotter.reset_camera._

- mock_main_window.view_3d_manager.plotter.reset_camera.assert_called_once()

### TestPluginInterface.test_enter_3d_viewer_mode
_enter_3d_viewer_mode delegates to ui_manager.enter_3d_viewer_mode._

- mock_main_window.ui_manager.enter_3d_viewer_mode.assert_called_once()

### TestPluginInterface.test_enter_3d_mode
_enter_3d_mode is an alias for enter_3d_viewer_mode._

- mock_main_window.ui_manager.enter_3d_viewer_mode.assert_called_once()

### TestPluginInterface.test_load_from_smiles
_load_from_smiles delegates to string_importer_manager and is safe when absent._

- mock_main_window.string_importer_manager.load_from_smiles.assert_called_once_with('C')

### TestPluginInterface.test_to_xyz_block
_to_xyz_block returns an XYZ string from current_mol, or None if no mol._

- assert xyz is not None
- assert 'C' in xyz
- assert '1.23' in xyz
- assert '4.56' in xyz
- assert '7.89' in xyz
- assert ctx.to_xyz_block() is None

### TestPluginInterface.test_set_bond_color_by_atoms
_set_bond_color_by_atoms looks up the bond and applies the color override._

- mock_mol.GetBondBetweenAtoms.assert_called_once_with(1, 2)
- mock_main_window.view_3d_manager.update_bond_color_override.assert_called_once_with(5, '#112233')
- mock_main_window.view_3d_manager.plotter.render.assert_called_once()

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
_register_toolbar_action stores the action in toolbar_actions._

- assert len(pm.toolbar_actions) == 1

### test_register_export_action
_register_export_action stores the action with its label in export_actions._

- assert len(pm.export_actions) == 1
- assert pm.export_actions[0]['label'] == 'Export as PDF'

### test_register_optimization_method
_register_optimization_method stores the callback under its method name._

- assert 'MMFF94' in pm.optimization_methods

### test_register_file_opener
_register_file_opener stores a callback under the given extension key._

- assert '.cif' in pm.file_openers

### test_register_file_opener_priority
_Higher priority opener should replace lower priority one in sorted order._

- assert pm.file_openers['.xyz'][0]['callback'] == cb_high

### test_register_analysis_tool
_register_analysis_tool appends a tool entry to analysis_tools._

- assert len(pm.analysis_tools) == 1

### test_register_save_handler
_register_save_handler stores the callback under the plugin name._

- assert 'TestPlugin' in pm.save_handlers

### test_register_load_handler
_register_load_handler stores the callback under the plugin name._

- assert 'TestPlugin' in pm.load_handlers

### test_register_3d_style
_register_3d_style stores the callback under the style name._

- assert 'Wireframe' in pm.custom_3d_styles

### test_register_document_reset_handler
_register_document_reset_handler appends to document_reset_handlers._

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

### test_install_and_discover_single_file
_install_plugin copies a .py file; discover_plugins returns its metadata._

- assert success
- assert (plugin_dir / 'test_plugin.py').exists()
- assert len(plugins) == 1
- assert p['name'] == 'Test Plugin'
- assert p['version'] == '1.0'
- assert p['description'] == 'A simple test plugin'
- assert p['status'] == 'Loaded'

### test_install_plugin_registers_menu_action
_A plugin that calls context.add_menu_action in initialize() should register it._

- assert len(pm.menu_actions) == 1
- assert pm.menu_actions[0]['plugin'] == 'Action Plugin'
- assert pm.menu_actions[0]['text'] == 'Test Action'

### test_install_zip_extracts_and_discovers
_install_plugin accepts a .zip file; extracted package is discovered._

- assert success
- assert (plugin_dir / 'MyPlugin' / '__init__.py').exists()
- assert any((p['name'] == 'Zipped Plugin' for p in plugins))

### test_compute_sha256
_Test SHA-256 calculation for files and directories._

- assert pm.compute_sha256('/non/existent/path') == 'N/A'
- assert sha1 != 'N/A'
- assert len(sha1) == 64
- assert sha_dir != 'N/A'
- assert len(sha_dir) == 64
- assert pm._sha256_for_file(str(f)) == 'N/A'
- assert pm._sha256_for_directory(str(d)) == 'N/A'

### TestPluginManagerExtended.test_imports_fallback
_PluginManager import does not crash when plugin_interface module is unavailable._


### TestPluginManagerExtended.test_get_set_main_window
_get_main_window returns the value set by set_main_window._

- assert pm.get_main_window() == 'mw'

### TestPluginManagerExtended.test_ensure_plugin_dir_error
_ensure_plugin_dir logs an error when os.makedirs raises OSError._

- mock_log.assert_called_with('Error creating plugin directory: test mkdir err')

### TestPluginManagerExtended.test_open_plugin_folder
_open_plugin_folder calls QDesktopServices.openUrl with the plugin folder URL._

- mock_open_url.assert_called_once()

### TestPluginManagerExtended.test_install_plugin_folder
_install_plugin handles a folder source, removing stale target when present._

- assert success
- assert 'package' in msg
- mock_rmtree.assert_called()
- mock_remove.assert_called()

### TestPluginManagerExtended.test_install_plugin_file_existing_dir
_install_plugin removes an existing directory before copying a file with the same name._

- mock_rmtree.assert_called()

### TestPluginManagerExtended.test_install_plugin_exception
_install_plugin returns (False, message) when an unexpected exception is raised._

- assert not success
- assert 'Install err' in msg

### TestPluginManagerExtended.test_zip_extraction
_install_plugin correctly extracts nested, flat, and edge-case ZIP structures._

- assert success
- assert success
- assert success
- assert success
- assert success

### TestPluginManagerExtended.test_discover_plugins_not_exists
_discover_plugins returns an empty list when the plugin directory does not exist._

- assert pm.discover_plugins() == []

### TestPluginManagerExtended.test_load_single_plugin_exceptions_and_stub
__load_single_plugin logs errors on RuntimeError and handles a None spec._

- assert 'Cat' in sys.modules
- assert 'Cat.SubCat' in sys.modules
- mock_log.assert_called()

### TestPluginManagerExtended.test_plugin_init_exceptions
_Plugins raising in initialize or autorun are recorded with an error status._

- assert plugin and 'Error (Init): INIT_ERR' in plugin['status']
- assert plugin2 and 'Error (Autorun): AUTO_ERR' in plugin2['status']

### TestPluginManagerExtended.test_load_plugin_version_tuple
__load_single_plugin converts a PLUGIN_VERSION tuple to a dotted string._

- assert pm.plugins[0]['version'] == '3.1.4'

### TestPluginManagerExtended.test_run_plugin_exceptions
_run_plugin shows a critical dialog when the plugin run() raises._

- mock_crit.assert_called_once()

### TestPluginManagerExtended.test_register_drop_handler
_register_drop_handler inserts handlers sorted by priority descending._

- assert pm.drop_handlers[0]['plugin'] == 'B'

### TestPluginManagerExtended.test_manager_api_helpers
_Status, undo, 3D-render, and camera-reset helpers delegate to main_window objects._

- mw.statusBar().showMessage.assert_called_with('test', 1000)
- mw.state_manager.push_undo_state.assert_called()
- mw.plotter.render.assert_called()
- mw.plotter.reset_camera.assert_called()
- mw.plotter.render.assert_called()

### TestPluginManagerExtended.test_get_selected_atom_indices_complex
_get_selected_atom_indices handles None main_window, missing atom_id, and type errors._

- assert pm.get_selected_atom_indices() == []
- assert indices == [0]
- assert pm.get_selected_atom_indices() == []
- assert pm.get_selected_atom_indices() == []

### TestPluginManagerExtended.test_register_get_window
_register_window stores and get_window retrieves a plugin's named window._

- assert pm.plugin_windows['P1']['w1'] == 'WIN'
- assert pm.get_window('P1', 'w1') == 'WIN'

### TestPluginManagerExtended.test_invoke_document_reset_handlers_logs_error
_invoke_document_reset_handlers logs an error when a handler raises._

- mock_log.assert_called()

### TestPluginManagerExtended.test_get_plugin_info_safe_exceptions_and_ast
_get_plugin_info_safe returns Unknown version for non-constant AST nodes and OSError._

- assert info['version'] == 'Unknown'

## tests/unit/test_plugin_manager_window.py

### test_init_and_refresh
_PluginManagerWindow initialises with correct title, row count, and status colours._

- assert window.windowTitle() == 'Plugin Manager'
- assert window.table.rowCount() == 3
- assert window.table.item(0, 0).text() == 'Loaded'
- assert window.table.item(0, 0).foreground().color() == Qt.GlobalColor.darkGreen
- assert window.table.item(1, 0).text() == 'Error'
- assert window.table.item(1, 0).foreground().color() == Qt.GlobalColor.red
- assert window.table.item(2, 0).text() == 'No Entry Point'
- assert window.table.item(2, 0).foreground().color() == Qt.GlobalColor.gray
- assert window.table.item(0, 4).text() == 'plugin1.py'

### test_refresh_relative_path_error
_When os.path.relpath raises ValueError, the filepath column falls back to basename._

- assert window.table.item(0, 4).text() == 'plugin3.py'

### test_update_button_state
_Remove button is disabled initially and enabled after a row is selected._

- assert not window.btn_remove.isEnabled()
- assert window.btn_remove.isEnabled()

### test_on_reload_main_window_present
_on_reload discovers plugins with the main window and shows an info message._

- mock_plugin_manager.discover_plugins.assert_called_with(mock_plugin_manager.main_window)
- mock_info.assert_called_once()
- mock_info.assert_not_called()

### test_on_reload_no_main_window
_on_reload calls discover_plugins without arguments when main_window is None._

- mock_plugin_manager.discover_plugins.assert_called_with()

### test_explore_plugins_online
_The Explore button opens the plugin explorer URL in the system browser._

- mock_open_url.assert_called_once()
- assert 'https://hiroyokoyama.github.io/moleditpy-plugins/explorer/' in mock_open_url.call_args[0][0].url()

### test_on_remove_plugin_no_selection
_on_remove_plugin shows a warning when no plugin row is selected._

- mock_warn.assert_called_with(window, 'Warning', 'Please select a plugin to remove.')

### test_on_remove_plugin_single_file
_on_remove_plugin removes a single-file plugin by calling os.remove._

- mock_remove.assert_called_with('/fake/plugins/plugin1.py')
- mock_info.assert_called()

### test_on_remove_plugin_package
_on_remove_plugin removes a package plugin by calling shutil.rmtree on its directory._

- mock_rmtree.assert_called_with('/fake/plugins/pkg_plugin')
- mock_info.assert_called()

### test_on_remove_plugin_error
_on_remove_plugin shows a critical error dialog when os.remove raises OSError._

- mock_critical.assert_called_with(window, 'Error', 'Failed to delete plugin: Remove error')

### test_on_remove_plugin_not_exists
_on_remove_plugin shows a warning when the plugin file does not exist._

- mock_warn.assert_called()

### test_show_plugin_details
_show_plugin_details displays an information dialog containing the plugin name._

- mock_info.assert_called_once()
- assert 'Test Plugin 1' in mock_info.call_args[0][2]

### test_drag_enter_event
_dragEnterEvent accepts URLs and ignores non-URL mime data._

- event.accept.assert_called_once()
- event.ignore.assert_called_once()

### test_drop_event_valid_files
_Dropping a .py file installs it and shows a success info message._

- mock_plugin_manager.install_plugin.assert_called_once_with('/some/file.py')
- mock_info.assert_called_once()
- assert 'Installed fine' in mock_info.call_args[0][2]

### test_drop_event_init_py_package
_Dropping an __init__.py file installs its parent folder as a package._

- mock_plugin_manager.install_plugin.assert_called_once_with('/some/folder')
- assert 'Error inst' in mock_info.call_args[0][2]

### test_drop_event_zip_file
_Dropping a .zip file installs it via install_plugin._

- mock_plugin_manager.install_plugin.assert_called_with('/some/file.zip')

### test_drop_event_pure_folder
_Dropping a directory installs it directly via install_plugin._

- mock_plugin_manager.install_plugin.assert_called_with('/some/plugin_folder')

## tests/unit/test_plugin_menu_manager.py

### TestConstruction.test_holds_init_manager_reference
_PluginMenuManager stores a reference to the init manager it received._

- assert mgr._im is im

### TestUpdatePluginMenu.test_does_nothing_when_no_plugin_manager
_update_plugin_menu does nothing when plugin_manager is None._

- menu.clear.assert_not_called()

### TestUpdatePluginMenu.test_clears_and_adds_manage_action
_update_plugin_menu clears the menu and adds at least the manage action._

- menu.clear.assert_called_once()
- assert menu.addAction.called or menu.addSeparator.called

### TestUpdatePluginMenu.test_discover_plugins_called
_update_plugin_menu calls discover_plugins on the plugin_manager._

- im.host.plugin_manager.discover_plugins.assert_called_once_with(im.host)

### TestUpdatePluginMenu.test_menu_action_appears_after_update
_A registered menu action is present in the menu after update_plugin_menu._

- im.host.menuBar.return_value.addMenu.assert_called_once_with('MyPlugin')
- added_menu.addAction.assert_called_once()

### TestRebuildPluginMenus.test_toolbar_and_export_populated_after_rebuild
_rebuild_plugin_menus wires toolbar and export actions into UI._

- im.plugin_toolbar.addAction.assert_called_once()
- im.export_button.menu.return_value.addAction.assert_called_once()

### TestRebuildPluginMenus.test_resets_separator_flag
_rebuild_plugin_menus resets the plugin_menubar_separator_added flag._

- assert im.plugin_menubar_separator_added is False

### TestRebuildPluginMenus.test_one_step_exception_does_not_stop_others
_A failing step must not prevent subsequent steps from running._

- assert 'boom' in call_order
- assert 'toolbar' in call_order
- assert 'export' in call_order
- assert 'style' in call_order

### TestRebuildPluginMenus.test_menu_cleanup_removes_tagged_actions
_Tagged plugin actions are stripped before rebuild._

- top_menu.removeAction.assert_called_with(tagged_action)

### TestAddRegisteredPluginActions.test_no_menu_actions_does_nothing
_add_registered_plugin_actions does nothing when no menu actions are registered._

- im.host.menuBar.return_value.addMenu.assert_not_called()

### TestAddRegisteredPluginActions.test_creates_new_top_level_menu_with_separator
_add_registered_plugin_actions creates a new top-level menu and adds a separator._

- im.host.menuBar.return_value.addSeparator.assert_called_once()
- im.host.menuBar.return_value.addMenu.assert_called_once_with('MyPlugin')
- new_menu.addAction.assert_called_once()

### TestAddRegisteredPluginActions.test_reuses_existing_top_level_menu
_add_registered_plugin_actions reuses an existing top-level menu rather than creating a new one._

- im.host.menuBar.return_value.addMenu.assert_not_called()
- existing_menu.addAction.assert_called_once()

### TestAddRegisteredPluginActions.test_separator_added_only_once_for_multiple_plugins
_The menubar separator is added exactly once even with multiple new menus._

- assert im.host.menuBar.return_value.addSeparator.call_count == 1

### TestAddRegisteredPluginActions.test_shortcut_applied_when_present
_add_registered_plugin_actions applies a keyboard shortcut when one is specified._

- assert isinstance(added_action, QAction)

### TestAddPluginToolbarActions.test_no_toolbar_attribute_does_nothing
_add_plugin_toolbar_actions does not raise when plugin_toolbar is missing._


### TestAddPluginToolbarActions.test_hides_toolbar_when_no_actions
_add_plugin_toolbar_actions hides and clears the toolbar when no actions are registered._

- im.plugin_toolbar.hide.assert_called_once()
- im.plugin_toolbar.clear.assert_called_once()

### TestAddPluginToolbarActions.test_shows_toolbar_and_adds_actions
_add_plugin_toolbar_actions shows the toolbar and adds an action for each registered entry._

- im.plugin_toolbar.show.assert_called_once()
- im.plugin_toolbar.addAction.assert_called_once()

### TestAddPluginToolbarActions.test_icon_set_when_file_exists
_add_plugin_toolbar_actions creates a QAction with an icon when the icon file exists._

- assert isinstance(added, QAction)

### TestIntegratePluginExportActions.test_no_export_actions_does_nothing
_integrate_plugin_export_actions does nothing when no export actions are registered._

- im.export_button.menu.return_value.addSeparator.assert_not_called()

### TestIntegratePluginExportActions.test_adds_actions_to_export_button_menu
_integrate_plugin_export_actions adds a separator and action to the export button menu._

- export_menu.addSeparator.assert_called_once()
- assert export_menu.addAction.called

### TestIntegratePluginExportActions.test_adds_actions_to_both_export_button_and_file_menu
_integrate_plugin_export_actions adds actions to both the export button and the File/Export submenu._

- assert im.export_button.menu.return_value.addAction.called
- assert export_submenu.addAction.called

### TestIntegratePluginAnalysisTools.test_no_analysis_menu_does_nothing
_integrate_plugin_analysis_tools does not raise when no Analysis menu exists._


### TestIntegratePluginAnalysisTools.test_no_tools_skips_separator
_integrate_plugin_analysis_tools skips the separator when no tools are registered._

- analysis_menu.addSeparator.assert_not_called()

### TestIntegratePluginAnalysisTools.test_adds_tools_to_analysis_menu
_integrate_plugin_analysis_tools adds a separator and action to the Analysis menu._

- analysis_menu.addSeparator.assert_called_once()
- analysis_menu.addAction.assert_called_once()

### TestUpdateStyleMenuWithPlugins.test_no_style_button_does_nothing
_update_style_menu_with_plugins does not raise when style_button is None._


### TestUpdateStyleMenuWithPlugins.test_no_custom_styles_does_nothing
_update_style_menu_with_plugins does nothing when no custom styles are registered._

- im.style_button.menu.return_value.addAction.assert_not_called()

### TestUpdateStyleMenuWithPlugins.test_adds_custom_style_actions
_update_style_menu_with_plugins adds a separator and one action per custom style._

- assert style_menu.addAction.call_count == 2
- assert style_menu.addSeparator.called

### TestUpdateStyleMenuWithPlugins.test_does_not_duplicate_existing_style
_update_style_menu_with_plugins skips a style already present in the menu._

- style_menu.addAction.assert_not_called()

### TestIntegratePluginFileOpeners.test_no_file_openers_does_nothing
_integrate_plugin_file_openers does nothing when no file openers are registered._

- im.import_menu.addSeparator.assert_not_called()

### TestIntegratePluginFileOpeners.test_no_import_menu_does_nothing
_integrate_plugin_file_openers does not raise when import_menu is None._


### TestIntegratePluginFileOpeners.test_adds_opener_actions
_integrate_plugin_file_openers adds a separator and action for each file opener._

- im.import_menu.addSeparator.assert_called_once()
- im.import_menu.addAction.assert_called_once()

### TestClearAllPluginActions.test_clears_plugin_menu
__clear_all_plugin_actions clears the plugin menu._

- plugin_menu.clear.assert_called_once()

### TestClearAllPluginActions.test_removes_tagged_actions_from_all_menus
__clear_all_plugin_actions removes all plugin_managed tagged actions from all menus._

- sub_menu.removeAction.assert_called_with(tagged)

## tests/unit/test_project_io.py

### test_save_project_no_data
_Verify error message when trying to save an empty project._

- io.statusBar().showMessage.assert_called_with('Error: Nothing to save.')

### test_save_project_overwrite_json
_Verify overwriting an existing JSON project file._

- assert os.path.exists(project_file)
- io.statusBar().showMessage.assert_any_call(f'Project saved to {project_file}')
- assert data['format'] == 'PME Project'

### test_save_project_overwrite_raw
_Verify overwriting an existing raw (pickle) project file._

- assert os.path.exists(raw_file)
- assert data['atoms'] == 'mock'

### test_save_project_redirect_to_save_as
_Verify that 'save' redirects to 'save as' if the current file is not a project file._

- assert mock_save_as.called

### test_load_raw_data_success
_Verify successful loading of a raw project file._

- assert mock_set_state.called
- assert io.host.init_manager.current_file_path == raw_file

### test_load_json_data_invalid_format
_Verify error handling for invalid JSON format in project files._

- assert mock_warn.called

### test_open_project_file_dispatch
_Verify dispatching to correct load method based on file extension._

- assert mock_json.called
- assert mock_raw.called

### test_save_as_json_trigger
_Verify triggering of 'save as' for JSON format._

- assert os.path.exists(save_path)

### test_load_raw_data_error_paths
_Verify error handling during raw data loading (file not found, corrupt)._

- io.statusBar().showMessage.assert_called_with('File not found: non_existent.pmeraw')
- io.statusBar().showMessage.assert_called_with('Invalid project file format: Corrupt')

### test_open_project_file_unsaved_check
_Verify that 'open project' checks for unsaved changes before proceeding._

- assert not mock_open.called

### test_save_project_io_error
_Verify handling of I/O errors during save._

- io.statusBar().showMessage.assert_called_with('File I/O error: Permission denied')

### test_load_json_data_version_mismatch
_Verify warning when loading a project from a newer software version._

- assert mock_info.called
- assert 'version 2.0' in mock_info.call_args[0][2]

### test_project_save_load_full_cycle
_Verify a complete save-load cycle for a project._

- assert os.path.exists(project_file)
- assert mock_load_json.called
- assert '2d_structure' in saved_data
- assert len(atoms) == 1
- assert atoms[0]['symbol'] == 'C'
- assert atoms[0]['charge'] == 1

### test_save_project_default_filename
_Verify that 'save as' suggests a reasonable default filename._

- assert suggested_path == expected
- assert args[2] == 'untitled'

### test_save_project_extension_enforcement
_Verify that correct file extensions are enforced when saving._

- assert io.host.init_manager.current_file_path == expected_path
- assert os.path.exists(expected_path)

### test_save_project_success_state_update
_Verify that application state is updated correctly after a successful save._

- assert io.host.state_manager.has_unsaved_changes is False
- assert io.host.init_manager.current_file_path == save_path
- assert io.host.update_window_title.called
- assert io.host.state_manager._saved_state is not None

### test_save_raw_data_no_data
_Verify error message when trying to save empty project._

- io.statusBar().showMessage.assert_called_with('Error: Nothing to save.')

### test_save_raw_data_success
_Verify successful saving via file dialog._

- assert os.path.exists(save_path)
- assert io.host.init_manager.current_file_path == save_path
- assert io.host.state_manager.has_unsaved_changes is False
- io.statusBar().showMessage.assert_called_with(f'Project saved to {save_path}')
- assert data == {'atoms': 'mock'}

### test_save_raw_data_cancel
_Verify that nothing happens if the user cancels the save dialog._

- io.statusBar().showMessage.assert_not_called()

### test_load_raw_data_dialog_success
_Verify loading via file dialog._

- io.host.state_manager.set_state_from_data.assert_called_with(sample_data)
- assert io.host.init_manager.current_file_path == load_path
- assert io.host.state_manager.has_unsaved_changes is False
- io.statusBar().showMessage.assert_called_with(f'Project loaded from {load_path}')

### test_load_raw_data_cancel
_Verify that nothing happens if the user cancels the load dialog._

- io.host.state_manager.set_state_from_data.assert_not_called()

### test_load_raw_data_io_error
_Verify handling of I/O errors during load._

- io.statusBar().showMessage.assert_called()
- assert 'Invalid project file format' in msg

## tests/unit/test_properties.py

### test_analysis_window_regular_mol
_Verify AnalysisWindow calculates and displays properties for regular molecules._

- assert formula_found
- assert 'C2H6O' in value

### test_analysis_window_xyz_derived
_Verify AnalysisWindow uses manual logic for XYZ-derived structures._

- assert 'C2HO' in formula_val
- assert not smiles_present

## tests/unit/test_radical_toggle.py

### test_radical_toggle_selected
_Test toggling radical on multiple selected atoms._

- assert a1.radical == 1
- assert a2.radical == 2
- assert scene.data.atoms[1]['radical'] == 1
- assert scene.data.atoms[2]['radical'] == 2
- a1.update_style.assert_called_once()
- a2.update_style.assert_called_once()
- scene.window.edit_actions_manager.push_undo_state.assert_called_once()
- event.accept.assert_called_once()

### test_radical_toggle_at_cursor
_Test toggling radical on an atom at the cursor when nothing is selected._

- assert a1.radical == 0
- assert scene.data.atoms[1]['radical'] == 0
- scene.window.edit_actions_manager.push_undo_state.assert_called_once()

### test_radical_toggle_no_target
_Test that nothing happens if no atoms are selected or at the cursor._

- scene.window.edit_actions_manager.push_undo_state.assert_not_called()
- event.accept.assert_not_called()

## tests/unit/test_ring_priority.py

### test_ring_priority_smaller_wins
_Verify that smaller rings are prioritized for double bond shift logic._

- assert shared_bond_item.is_in_ring is True
- assert shared_bond_item.ring_center == five_center

## tests/unit/test_scene_advanced.py

### test_right_click_bond_deletion
_Test standard right-click deletion on a bond._

- assert len(data.bonds) == 1
- assert len(data.bonds) == 0
- assert a1_id in data.atoms
- assert a2_id in data.atoms
- assert bond_item.scene() is None

### test_drag_and_drop_atom
_Test moving an atom via drag-and-drop._

- assert a1_item.pos() == new_pos
- assert data.atoms[a1_id]['pos'] == (new_pos.x(), new_pos.y())
- assert len(window.edit_actions_manager.undo_stack) > 0

### test_delete_mixed_selection
_Test deleting a selection containing both atoms and bonds._

- assert a1_id not in data.atoms
- assert len(data.bonds) == 0
- assert a2_id in data.atoms

### test_undo_redo
_Test undo/redo integration via scene modifications._

- assert len(window.edit_actions_manager.undo_stack) == 0
- assert len(window.edit_actions_manager.undo_stack) == initial_len + 1

### test_atom_fusing_enabled
_Test that drawing a bond near an existing atom snaps/fuses to it when enabled._

- assert len(data.atoms) == 2
- assert len(data.bonds) == 1

### test_atom_fusing_disabled
_Test that drawing a bond near an existing atom does not snap/fuse when disabled._

- assert len(data.atoms) == 3
- assert len(data.bonds) == 1

### test_atom_mode_short_click_near_atom_does_not_fuse
_Atom-mode short clicks near an atom should create a new atom, not edit the nearby one._

- assert data.atoms[existing_id]['symbol'] == 'C'
- assert len(data.atoms) == 2

### test_benzene_terminal_120_deg_alignment
_Test that pressing 4 at a terminal atom aligns the benzene ring at 120 degrees._

- assert len(data.atoms) == 7
- assert found_expected_pos

### test_click_selected_atom_deselects
_Test that clicking (without dragging) on a selected atom deselects it._

- assert scene.was_selected_on_press is True
- assert a1_item.isSelected() is False

### test_drag_selected_atom_keeps_selected
_Test that dragging a selected atom moves it but keeps it selected._

- assert a1_item.pos() == new_pos
- assert a1_item.isSelected() is True

## tests/unit/test_scene_extended.py

### test_scene_keypress_modes
_Letter key presses switch the scene to the corresponding atom placement mode._

- assert scene.window.ui_manager.activate_select_mode.called
- assert scene.mode == expected_mode

### test_scene_keypress_special_symbols
_Shift+letter shortcuts set two-letter atom modes like Cl, Br, Si._

- assert scene.mode == 'atom_Cl'
- assert scene.mode == 'atom_Br'
- assert scene.mode == 'atom_Si'

### test_scene_keypress_delete
_Delete key calls delete_items on the current selection._

- mock_delete.assert_called()

### test_scene_maintenance_methods
_Scene maintenance methods run without error and clear expected state._

- assert atom_item.has_problem is False

### test_scene_queries
_find_atom_near returns the atom closest to the query point._

- assert scene.find_atom_near(QPointF(50, 0)) == a2

### test_scene_update_connected_bonds
_update_connected_bonds calls update_position on each connected bond._

- assert a.bonds[0].update_position.called

### test_scene_leave_event
_Leaving the scene hides the template preview._

- assert scene.template_preview.hide.called

### test_scene_update_bond_stereo
_update_bond_stereo updates both the data dict and the bond item._

- assert scene.data.bonds[10, 20]['stereo'] == 1
- assert bond.set_stereo.called

### test_scene_mouse_drag_create_bond_existing_atoms
_Drag from one existing atom to another creates a bond between them._

- assert (aid1, aid2) in scene.data.bonds

### test_scene_mouse_click_create_single_atom
_Click on blank space in atom mode creates a new atom and pushes undo._

- assert any((a['symbol'] == 'O' for a in scene.data.atoms.values()))
- assert mock_parser_host.edit_actions_manager.push_undo_state.called

### test_scene_right_click_bond_delete
_Test right-click deletion on a bond._

- mock_del.assert_called()
- assert len(mock_parser_host.data.bonds) == 0
- assert mock_parser_host.edit_actions_manager.push_undo_state.called

### test_scene_drag_and_drop_atom
_Test moving an atom via drag-and-drop._

- assert atom_item.pos() == new_pos
- assert mock_parser_host.data.atoms[aid]['pos'] == (new_pos.x(), new_pos.y())
- assert mock_parser_host.edit_actions_manager.push_undo_state.called

### test_scene_delete_mixed_selection
_Test deleting a selection containing both atoms and bonds._

- assert aid1 not in mock_parser_host.data.atoms
- assert len(mock_parser_host.data.bonds) == 0
- assert aid2 in mock_parser_host.data.atoms

## tests/unit/test_scene_extended_logic.py

### test_scene_ez_toggle_logic
_Test E/Z stereo toggling logic in mouseReleaseEvent._

- assert bond.stereo == 3
- assert data.bonds[0, 1]['stereo'] == 3

### test_scene_item_deletion_path
_Test item deletion logic in mouseReleaseEvent._

- assert scene.delete_items.called

### test_scene_bond_inversion
_Test bond inversion logic in mouseReleaseEvent._

- assert data.remove_bond.called
- assert (1, 0) in data.bonds

### test_update_user_template_preview
_Test Template preview logic._

- assert 'points' in scene.template_context

### test_benzene_template_rotation_logic
_Test benzene template rotation alignment._

- assert data.add_atom.called

## tests/unit/test_scene_interactions.py

### test_scene_toggle_radical
_Radical mode cycles atom radical count through 0→1→2→0 on successive clicks._

- assert atom_item.radical == 1
- assert atom_item.radical == 2
- assert atom_item.radical == 0

### test_scene_toggle_charge
_Charge plus/minus modes increment and decrement atom charge correctly._

- assert atom_item.charge == 1
- assert atom_item.charge == 0

### test_add_benzene_fragment
_add_molecule_fragment creates 6 atoms and 6 bonds with correct alternating order._

- assert len(scene.data.atoms) == 6
- assert len(scene.data.bonds) == 6
- assert len([b for b in scene.data.bonds.values() if b['order'] == 2]) == 3

### test_benzene_fusion_rotation
_Benzene fusion reuses existing edge atoms and preserves bond order._

- assert len(scene.data.atoms) == 6
- assert eb.order == 2

### test_delete_selected_items
_Right-click on a selected atom triggers delete_items._

- mock_delete.assert_called()

### test_double_click_select_component
_Double-click selects all atoms in the connected component, not isolated ones._

- assert scene.atom_items[id1].setSelected.called
- assert scene.atom_items[id2].setSelected.called
- assert scene.atom_items[id3].setSelected.called
- assert not atom_iso.setSelected.called

### test_scene_key_event_dispatch
_Test key events like '4' (template), '.' (radical), +/- (charge)._

- mock_parser_host.ui_manager.set_mode_and_update_toolbar.assert_called_with('template_benzene')
- assert atom_item.charge == 1
- assert scene.data.atoms[aid]['charge'] == 1

### test_scene_update_template_preview_logic
_Test update_template_preview with different targets._

- scene.template_preview.set_geometry.assert_called()
- assert len(points) == 6
- scene.template_preview.set_geometry.assert_called()
- assert 'items' in scene.template_context
- assert scene.template_context['items'][0] == atom_item

### test_scene_drag_create_bond_sequence
_Test the full mouse press -> move -> release sequence for creating a bond._

- assert bond.order == 1
- assert getattr(scene, 'start_atom', None) == a1

### test_scene_bond_snapping_distance
_Verify that mouse interactions respect the configured bond_snapping_distance_2d setting._

- mock_find_near.assert_called_with(QPointF(10, 10), tol=25.0)

## tests/unit/test_settings_2d_tab.py

### test_init_uses_default_colors
_Settings2DTab initialises color fields from DEFAULT_SETTINGS._

- assert tab.current_bg_color_2d == DEFAULT_SETTINGS['background_color_2d']
- assert tab.current_bond_color_2d == DEFAULT_SETTINGS['bond_color_2d']

### test_update_ui_sets_colors
_update_ui() reflects background and bond color changes in internal state._

- assert tab.current_bg_color_2d == '#112233'
- assert tab.current_bond_color_2d == '#aabbcc'

### test_update_ui_sets_sliders
_update_ui() updates all bond and font slider values correctly._

- assert tab.bond_width_2d_slider.value() == 30
- assert tab.bond_spacing_double_2d_slider.value() == 40
- assert tab.bond_spacing_triple_2d_slider.value() == 50
- assert tab.bond_wedge_width_2d_slider.value() == 80
- assert tab.bond_dash_count_2d_slider.value() == 12
- assert tab.atom_font_size_2d_slider.value() == 24
- assert tab.template_fusing_distance_2d_slider.value() == 18
- assert tab.template_snapping_distance_2d_slider.value() == 20
- assert tab.bond_snapping_distance_2d_slider.value() == 16

### test_update_ui_sets_cap_style
_update_ui() reflects bond_cap_style_2d in the combo box._

- assert tab.bond_cap_style_2d_combo.currentText() == 'Flat'

### test_update_ui_sets_use_bond_color_checkbox
_update_ui() reflects atom_use_bond_color_2d in the checkbox state._

- assert tab.atom_use_bond_color_2d_checkbox.isChecked() is True

### test_get_settings_roundtrip
_get_settings() returns values matching what was loaded via update_ui()._

- assert result['background_color_2d'] == DEFAULT_SETTINGS['background_color_2d']
- assert result['bond_color_2d'] == DEFAULT_SETTINGS['bond_color_2d']
- assert abs(result['bond_width_2d'] - DEFAULT_SETTINGS['bond_width_2d']) < 0.05
- assert result['bond_cap_style_2d'] == DEFAULT_SETTINGS['bond_cap_style_2d']
- assert result['bond_dash_count_2d'] == DEFAULT_SETTINGS['bond_dash_count_2d']
- assert result['atom_font_size_2d'] == DEFAULT_SETTINGS['atom_font_size_2d']
- assert result['atom_use_bond_color_2d'] == DEFAULT_SETTINGS['atom_use_bond_color_2d']
- assert result['template_fusing_enabled_2d'] == DEFAULT_SETTINGS['template_fusing_enabled_2d']
- assert abs(result['template_fusing_distance_2d'] - DEFAULT_SETTINGS['template_fusing_distance_2d']) < 0.05
- assert abs(result['template_snapping_distance_2d'] - DEFAULT_SETTINGS['template_snapping_distance_2d']) < 0.05
- assert abs(result['bond_snapping_distance_2d'] - DEFAULT_SETTINGS['bond_snapping_distance_2d']) < 0.05

### test_get_settings_returns_all_keys
_get_settings() returns a dict containing all expected 2D settings keys._

- assert expected_keys == set(result.keys())

### test_atom_font_style_defaults
_Bold/Italic/Underline buttons reflect DEFAULT_SETTINGS on init._

- assert tab.atom_font_bold_2d_btn.isChecked() == DEFAULT_SETTINGS['atom_font_bold_2d']
- assert tab.atom_font_italic_2d_btn.isChecked() == DEFAULT_SETTINGS['atom_font_italic_2d']
- assert tab.atom_font_underline_2d_btn.isChecked() == DEFAULT_SETTINGS['atom_font_underline_2d']

### test_atom_font_style_update_ui
_update_ui() sets Bold/Italic/Underline button state from settings dict._

- assert tab.atom_font_bold_2d_btn.isChecked() is False
- assert tab.atom_font_italic_2d_btn.isChecked() is True
- assert tab.atom_font_underline_2d_btn.isChecked() is True

### test_atom_font_style_get_settings_roundtrip
_get_settings() returns Bold/Italic/Underline matching what was set via update_ui()._

- assert result['atom_font_bold_2d'] is False
- assert result['atom_font_italic_2d'] is True
- assert result['atom_font_underline_2d'] is True

### test_pick_bg_color_updates_on_valid
_Selecting a valid color from the dialog updates current_bg_color_2d._

- assert tab.current_bg_color_2d == '#ff0000'

### test_pick_bg_color_no_change_on_invalid
_Cancelling the color dialog leaves current_bg_color_2d unchanged._

- assert tab.current_bg_color_2d == original

### test_pick_bond_color_updates_on_valid
_Selecting a valid color from the dialog updates current_bond_color_2d._

- assert tab.current_bond_color_2d == '#00ff00'

### test_reset_to_defaults
_reset_to_defaults() restores all settings to DEFAULT_SETTINGS values._

- assert tab.current_bg_color_2d == '#123456'
- assert tab.template_fusing_enabled_2d_checkbox.isChecked() is False
- assert tab.current_bg_color_2d == DEFAULT_SETTINGS['background_color_2d']
- assert tab.bond_cap_style_2d_combo.currentText() == DEFAULT_SETTINGS['bond_cap_style_2d']
- assert tab.template_fusing_enabled_2d_checkbox.isChecked() == DEFAULT_SETTINGS['template_fusing_enabled_2d']
- assert tab.template_fusing_distance_2d_slider.value() == int(DEFAULT_SETTINGS['template_fusing_distance_2d'])
- assert tab.template_snapping_distance_2d_slider.value() == int(DEFAULT_SETTINGS['template_snapping_distance_2d'])
- assert tab.bond_snapping_distance_2d_slider.value() == int(DEFAULT_SETTINGS['bond_snapping_distance_2d'])

### test_template_fusing_checkbox_disables_slider
_Template fusing checkbox enables/disables the distance slider and label._

- assert tab.template_fusing_distance_2d_slider.isEnabled() is True
- assert tab.template_fusing_distance_2d_label.isEnabled() is True
- assert tab.template_fusing_distance_2d_slider.isEnabled() is False
- assert tab.template_fusing_distance_2d_label.isEnabled() is False

## tests/unit/test_settings_3d_tabs.py

### test_scene_tab_init
_Settings3DSceneTab initialises with the default background color._

- assert tab.current_bg_color == DEFAULT_SETTINGS['background_color']

### test_scene_tab_update_ui
_update_ui applies all provided settings to the scene tab widgets._

- assert tab.current_bg_color == '#112233'
- assert tab.axes_checkbox.isChecked() is False
- assert tab.light_checkbox.isChecked() is False
- assert tab.intensity_slider.value() == 50
- assert tab.specular_slider.value() == 10
- assert tab.spec_power_slider.value() == 40
- assert tab.projection_combo.currentText() == 'Orthographic'

### test_scene_tab_get_settings_keys
_get_settings returns a dict with all expected scene setting keys._

- assert expected_keys == set(result.keys())

### test_scene_tab_roundtrip
_Settings round-trip: update_ui then get_settings recovers the original values._

- assert result['background_color'] == DEFAULT_SETTINGS['background_color']
- assert result['show_3d_axes'] == DEFAULT_SETTINGS['show_3d_axes']
- assert result['projection_mode'] == DEFAULT_SETTINGS['projection_mode']
- assert abs(result['light_intensity'] - DEFAULT_SETTINGS['light_intensity']) < 0.01
- assert abs(result['specular'] - DEFAULT_SETTINGS['specular']) < 0.01
- assert result['specular_power'] == DEFAULT_SETTINGS['specular_power']

### test_scene_tab_pick_color_valid
__select_color updates current_bg_color when a valid color is chosen._

- assert tab.current_bg_color == '#abcdef'

### test_scene_tab_pick_color_invalid_no_change
__select_color leaves current_bg_color unchanged when an invalid color is returned._

- assert tab.current_bg_color == original

### test_scene_tab_reset_to_defaults
_reset_to_defaults restores projection_mode to the default setting._

- assert tab.projection_combo.currentText() == 'Orthographic'
- assert tab.projection_combo.currentText() == DEFAULT_SETTINGS['projection_mode']

### test_model_tab_ball_stick_has_atom_scale
_ball_stick SettingsModelTab creates an atom_scale_slider widget._

- assert hasattr(tab, 'atom_scale_slider')

### test_model_tab_ball_stick_has_bond_radius
_ball_stick SettingsModelTab creates a bond_radius_slider widget._

- assert hasattr(tab, 'bond_radius_slider')

### test_model_tab_ball_stick_has_color_options
_ball_stick SettingsModelTab creates bond_color_button and use_cpk_checkbox._

- assert hasattr(tab, 'bond_color_button')
- assert hasattr(tab, 'use_cpk_checkbox')

### test_model_tab_ball_stick_update_ui
_update_ui sets all ball-and-stick sliders and checkbox from settings._

- assert tab.atom_scale_slider.value() == 150
- assert tab.bond_radius_slider.value() == 20
- assert tab.res_slider.value() == 20
- assert tab.use_cpk_checkbox.isChecked() is True

### test_model_tab_ball_stick_get_settings_keys
_get_settings returns all expected ball-and-stick setting keys._

- assert expected_keys == set(result.keys())

### test_model_tab_ball_stick_roundtrip
_ball_stick settings round-trip recovers atom scale, bond radius, and resolution._

- assert abs(result['ball_stick_atom_scale'] - DEFAULT_SETTINGS['ball_stick_atom_scale']) < 0.01
- assert abs(result['ball_stick_bond_radius'] - DEFAULT_SETTINGS['ball_stick_bond_radius']) < 0.01
- assert result['ball_stick_resolution'] == DEFAULT_SETTINGS['ball_stick_resolution']

### test_model_tab_cpk_has_atom_scale_no_bond_radius
_cpk SettingsModelTab has atom_scale_slider but not bond_radius_slider._

- assert hasattr(tab, 'atom_scale_slider')
- assert not hasattr(tab, 'bond_radius_slider')

### test_model_tab_cpk_get_settings_keys
_cpk get_settings returns cpk_atom_scale and cpk_resolution, not ball-stick keys._

- assert 'cpk_atom_scale' in result
- assert 'cpk_resolution' in result
- assert 'ball_stick_bond_color' not in result

### test_model_tab_cpk_roundtrip
_cpk settings round-trip recovers atom scale and resolution._

- assert abs(result['cpk_atom_scale'] - DEFAULT_SETTINGS['cpk_atom_scale']) < 0.01
- assert result['cpk_resolution'] == DEFAULT_SETTINGS['cpk_resolution']

### test_model_tab_wireframe_no_atom_scale
_wireframe SettingsModelTab has bond_radius_slider but not atom_scale_slider._

- assert not hasattr(tab, 'atom_scale_slider')
- assert hasattr(tab, 'bond_radius_slider')

### test_model_tab_wireframe_get_settings_keys
_wireframe get_settings returns wireframe_bond_radius, resolution, and offset keys._

- assert 'wireframe_bond_radius' in result
- assert 'wireframe_resolution' in result
- assert 'wireframe_double_bond_offset_factor' in result

### test_model_tab_wireframe_roundtrip
_wireframe settings round-trip recovers bond radius and resolution._

- assert abs(result['wireframe_bond_radius'] - DEFAULT_SETTINGS['wireframe_bond_radius']) < 0.01
- assert result['wireframe_resolution'] == DEFAULT_SETTINGS['wireframe_resolution']

### test_model_tab_stick_get_settings_keys
_stick get_settings returns stick_bond_radius, resolution, and offset keys._

- assert 'stick_bond_radius' in result
- assert 'stick_resolution' in result
- assert 'stick_double_bond_offset_factor' in result

### test_model_tab_stick_roundtrip
_stick settings round-trip recovers bond radius and resolution._

- assert abs(result['stick_bond_radius'] - DEFAULT_SETTINGS['stick_bond_radius']) < 0.01
- assert result['stick_resolution'] == DEFAULT_SETTINGS['stick_resolution']

### test_model_tab_ball_stick_pick_bond_color
__pick_bond_color updates current_bond_color when a valid color is chosen._

- assert tab.current_bond_color == '#ff00ff'

## tests/unit/test_settings_dialog.py

### test_init_creates_seven_tabs
_SettingsDialog initialises with exactly seven tabs._

- assert dialog.tab_widget.count() == 7

### test_init_tab_labels
_Tab titles include the expected section names._

- assert '2D Settings' in tab_titles
- assert '3D Scene' in tab_titles
- assert 'Ball & Stick' in tab_titles
- assert 'Other' in tab_titles

### test_init_window_title
_SettingsDialog window title is 'Settings'._

- assert dialog.windowTitle() == 'Settings'

### test_update_ui_from_settings_sets_all_tabs
_update_ui_from_settings propagates settings values to each tab widget._

- assert dialog.tab_2d.current_bg_color_2d == '#aabbcc'

### test_get_settings_aggregates_all_tabs
_get_settings aggregates keys from every settings tab._

- assert 'background_color_2d' in result
- assert 'background_color' in result
- assert 'ball_stick_atom_scale' in result
- assert 'cpk_atom_scale' in result
- assert 'wireframe_bond_radius' in result
- assert 'stick_bond_radius' in result
- assert 'skip_chemistry_checks' in result

### test_reset_current_tab_calls_reset_on_active_tab
_reset_current_tab calls reset_to_defaults on the currently visible tab._

- mock_reset.assert_called_once()

### test_reset_current_tab_shows_info_message
_reset_current_tab shows an information message naming the reset tab._

- mock_info.assert_called_once()
- assert '2D Settings' in args[2]

### test_reset_all_settings_yes_resets_and_applies
_reset_all_settings resets UI and applies when user confirms Yes._

- mock_update.assert_called_once_with(dialog.default_settings)
- mock_apply.assert_called_once()

### test_reset_all_settings_no_does_nothing
_reset_all_settings is a no-op when user answers No._

- mock_update.assert_not_called()

### test_apply_settings_no_parent_returns_early
_apply_settings is a no-op when parent_window is None._


### test_apply_settings_updates_parent_settings
_apply_settings pushes current tab values into the parent's settings dict._

- assert dialog.parent_window.init_manager.settings['background_color_2d'] == '#123456'

### test_apply_settings_fires_all_side_effects
_apply_settings calls save_settings, apply_3d_settings, update_cpk_colors, and status bar._

- p.init_manager.save_settings.assert_called()
- p.view_3d_manager.apply_3d_settings.assert_called()
- p.init_manager.update_cpk_colors_from_settings.assert_called()
- p.statusBar.return_value.showMessage.assert_called_with('Settings applied successfully')

### test_apply_settings_redraws_molecule_when_present
_apply_settings calls draw_molecule_3d when a molecule is currently loaded._

- parent.view_3d_manager.draw_molecule_3d.assert_called_with(mol)

### test_apply_settings_skips_redraw_when_no_molecule
_apply_settings skips draw_molecule_3d when no molecule is loaded._

- parent.view_3d_manager.draw_molecule_3d.assert_not_called()

### test_apply_settings_updates_2d_scene_background
_apply_settings calls setBackgroundBrush on the 2D scene._

- scene.setBackgroundBrush.assert_called()

### test_accept_applies_then_closes
_accept() calls apply_settings before delegating to QDialog.accept._

- mock_apply.assert_called_once()
- mock_super_accept.assert_called_once()

### test_get_settings_roundtrip_defaults
_get_settings after loading DEFAULT_SETTINGS round-trips the default values._

- assert abs(result['ball_stick_atom_scale'] - DEFAULT_SETTINGS['ball_stick_atom_scale']) < 0.01
- assert result['background_color'] == DEFAULT_SETTINGS['background_color']
- assert result['skip_chemistry_checks'] == DEFAULT_SETTINGS['skip_chemistry_checks']

## tests/unit/test_settings_other_tab.py

### test_init_creates_all_controls
_SettingsOtherTab creates all expected control widgets on init._

- assert hasattr(tab, 'skip_chem_checks_checkbox')
- assert hasattr(tab, 'always_ask_charge_checkbox')
- assert hasattr(tab, 'kekule_3d_checkbox')
- assert hasattr(tab, 'aromatic_circle_checkbox')
- assert hasattr(tab, 'aromatic_torus_thickness_slider')

### test_update_ui_sets_checkboxes
_update_ui sets checkbox states from provided settings dict._

- assert tab.skip_chem_checks_checkbox.isChecked() is True
- assert tab.always_ask_charge_checkbox.isChecked() is True
- assert tab.kekule_3d_checkbox.isChecked() is False
- assert tab.aromatic_circle_checkbox.isChecked() is False

### test_update_ui_sets_torus_thickness
_update_ui maps aromatic_torus_thickness_factor to slider integer value._

- assert tab.aromatic_torus_thickness_slider.value() == 120

### test_get_settings_returns_correct_keys
_get_settings returns a dict with exactly the expected setting keys._

- assert expected_keys == set(result.keys())

### test_get_settings_roundtrip_defaults
_get_settings after update_ui(DEFAULT_SETTINGS) round-trips the default values._

- assert result['skip_chemistry_checks'] == DEFAULT_SETTINGS['skip_chemistry_checks']
- assert result['display_kekule_3d'] == DEFAULT_SETTINGS['display_kekule_3d']
- assert result['display_aromatic_circles_3d'] == DEFAULT_SETTINGS['display_aromatic_circles_3d']
- assert abs(result['aromatic_torus_thickness_factor'] - DEFAULT_SETTINGS['aromatic_torus_thickness_factor']) < 0.01

### test_kekule_toggled_disables_aromatic_circle
_Checking Kekulé disables the aromatic circle checkbox._

- assert tab.aromatic_circle_checkbox.isEnabled() is False

### test_kekule_untoggled_enables_aromatic_circle
_Unchecking Kekulé re-enables the aromatic circle checkbox._

- assert tab.aromatic_circle_checkbox.isEnabled() is True

### test_aromatic_toggled_disables_kekule
_Checking aromatic circles disables the Kekulé checkbox._

- assert tab.kekule_3d_checkbox.isEnabled() is False

### test_aromatic_untoggled_enables_kekule
_Unchecking aromatic circles re-enables the Kekulé checkbox._

- assert tab.kekule_3d_checkbox.isEnabled() is True

### test_update_ui_kekule_true_disables_aromatic
_update_ui with display_kekule_3d=True disables the aromatic circle checkbox._

- assert tab.kekule_3d_checkbox.isChecked() is True
- assert tab.aromatic_circle_checkbox.isEnabled() is False

### test_update_ui_aromatic_true_disables_kekule
_update_ui with display_aromatic_circles_3d=True disables the Kekulé checkbox._

- assert tab.aromatic_circle_checkbox.isChecked() is True
- assert tab.kekule_3d_checkbox.isEnabled() is False

### test_reset_to_defaults
_reset_to_defaults restores settings to default values._

- assert result['skip_chemistry_checks'] == DEFAULT_SETTINGS['skip_chemistry_checks']
- assert abs(result['aromatic_torus_thickness_factor'] - DEFAULT_SETTINGS['aromatic_torus_thickness_factor']) < 0.01

### test_sync_slider_from_spinbox
_Setting the spinbox value syncs the torus thickness slider accordingly._

- assert tab.aromatic_torus_thickness_slider.value() == 150

### test_log_to_file_default_unchecked
_log_to_file checkbox is unchecked by default._

- assert tab.log_to_file_checkbox.isChecked() is False
- assert tab.log_level_debug_checkbox.isChecked() is False

### test_log_settings_roundtrip
_Enabling log_to_file and log_level_debug round-trips through get_settings._

- assert result['log_to_file'] is True
- assert result['log_level_debug'] is True

### test_log_settings_off_roundtrip
_Disabling both log settings returns False in get_settings._

- assert result['log_to_file'] is False
- assert result['log_level_debug'] is False

## tests/unit/test_settings_tab_base.py

### test_init_stores_default_settings
___init__ stores the provided default_settings dict for later use._

- assert tab.default_settings is DEFAULT_SETTINGS

### test_create_separator_returns_hline
__create_separator returns an HLine QFrame with Sunken shadow._

- assert isinstance(sep, QFrame)
- assert sep.frameShape() == QFrame.Shape.HLine
- assert sep.frameShadow() == QFrame.Shadow.Sunken

### test_create_slider_range
__create_slider returns a QSlider with the specified min/max range._

- assert isinstance(slider, QSlider)
- assert slider.minimum() == 10
- assert slider.maximum() == 200

### test_create_slider_float_spin_updates
_Slider and float spinbox stay in sync bidirectionally._

- assert spin.value() == 5.0
- assert slider.value() == 80

### test_create_slider_int_spin_updates
_Slider and integer spinbox stay in sync bidirectionally._

- assert spin.value() == 10
- assert slider.value() == 15

### test_wrap_layout_returns_layout_with_children
__wrap_layout returns an HBoxLayout containing the provided widgets._

- assert isinstance(layout, QHBoxLayout)
- assert layout.indexOf(slider) != -1
- assert layout.indexOf(label) != -1

### test_reset_to_defaults_calls_update_ui
_reset_to_defaults calls update_ui with the stored default_settings._

- tab.update_ui.assert_called_once_with(DEFAULT_SETTINGS)

### test_update_ui_not_implemented
_SettingsTabBase.update_ui raises NotImplementedError on unsubclassed base._


### test_get_settings_not_implemented
_SettingsTabBase.get_settings raises NotImplementedError on unsubclassed base._


## tests/unit/test_sip_isdeleted_safe.py

### TestSipIsDeletedSafe.test_none_is_treated_as_deleted
_None input is treated as a deleted object, returning True._

- assert sip_isdeleted_safe(None) is True

### TestSipIsDeletedSafe.test_returns_false_when_sip_unavailable
_Returns False (not deleted) when the sip checker is not available._

- assert sip_isdeleted_safe(object()) is False

### TestSipIsDeletedSafe.test_returns_true_when_sip_reports_deleted
_Returns True when sip confirms the object is deleted._

- assert sip_isdeleted_safe(object()) is True

### TestSipIsDeletedSafe.test_returns_false_when_sip_reports_not_deleted
_Returns False when sip confirms the object is still alive._

- assert sip_isdeleted_safe(object()) is False

### TestSipIsDeletedSafe.test_returns_false_on_runtime_error
_Returns False conservatively when sip raises RuntimeError._

- assert sip_isdeleted_safe(object()) is False

### TestSipIsDeletedSafe.test_returns_false_on_attribute_error
_Returns False conservatively when sip raises AttributeError._

- assert sip_isdeleted_safe(object()) is False

### TestSipIsDeletedSafe.test_returns_false_on_type_error
_Returns False conservatively when sip raises TypeError._

- assert sip_isdeleted_safe(object()) is False

### TestSipIsDeletedSafe.test_fallback_to_standard_sip
_Falls back to standard sip.isdeleted when PyQt6.sip is unavailable._

- assert sds.sip_isdeleted_safe(object()) is True

### TestSipIsDeletedSafe.test_fallback_when_both_sip_packages_missing
_Returns False and sets _sip_isdeleted=None when both sip packages are absent._

- assert sds.sip_isdeleted_safe(object()) is False
- assert sds._sip_isdeleted is None

## tests/unit/test_slider_logic.py

### test_angle_dialog_wrapping
_Angle input > 180 wraps back into the valid signed range._

- assert dialog.angle_input.text() == '-170.00'
- dialog.apply_geometry_update.assert_called_once_with(-170.0)

### test_dihedral_dialog_wrapping
_Dihedral input outside [-180, 180] wraps back into the valid signed range._

- assert dialog.dihedral_input.text() == '160.00'
- dialog.apply_geometry_update.assert_called_once_with(160.0)

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

- mock_parser_host.statusBar().showMessage.assert_called()
- assert last_msg.startswith('Invalid SMILES:')

### test_smiles_empty_shows_error
_Empty SMILES should trigger error message._

- mock_parser_host.statusBar().showMessage.assert_called()

### test_inchi_ethanol
_InChI import should match RDKit reference atom/bond count._

- assert len(data.atoms) == ref.GetNumAtoms()
- assert len(data.bonds) == ref.GetNumBonds()

### test_inchi_caffeine_formula
_Caffeine InChI should produce correct molecular formula._

- assert rdMolDescriptors.CalcMolFormula(mol) == rdMolDescriptors.CalcMolFormula(ref)

### test_inchi_invalid_shows_error
_Invalid InChI should show error, not crash._

- mock_parser_host.statusBar().showMessage.assert_called()
- assert last_msg.startswith('Invalid InChI:')

## tests/unit/test_system_utils.py

### TestWindowsTheme.test_windows_dark_when_val_zero
_AppsUseLightTheme=0 maps to dark theme on Windows._

- assert detect_system_theme() == 'dark'

### TestWindowsTheme.test_windows_light_when_val_one
_AppsUseLightTheme=1 maps to light theme on Windows._

- assert detect_system_theme() == 'light'

### TestWindowsTheme.test_windows_none_when_winreg_unavailable
_If winreg is None (non-Windows build), Windows branch is skipped._

- assert detect_system_theme() is None

### TestWindowsTheme.test_windows_oserror_returns_none
_OSError in winreg must be suppressed and return None._

- assert detect_system_theme() is None

### TestMacOSTheme.test_macos_always_light
_macOS always returns light theme._

- assert detect_system_theme() == 'light'

### TestLinuxTheme.test_linux_gnome_color_scheme_dark
_GNOME color-scheme prefer-dark maps to dark theme on Linux._

- assert detect_system_theme() == 'dark'

### TestLinuxTheme.test_linux_gnome_color_scheme_light
_GNOME color-scheme prefer-light maps to light theme on Linux._

- assert detect_system_theme() == 'light'

### TestLinuxTheme.test_linux_gtk_theme_dark_fallback
_When color-scheme returns nothing useful, fall back to gtk-theme name._

- assert detect_system_theme() == 'dark'

### TestLinuxTheme.test_linux_gtk_theme_no_dark_keyword
_gtk-theme without '-dark' should not return dark._

- assert detect_system_theme() is None

### TestLinuxTheme.test_linux_gsettings_not_found_returns_none
_FileNotFoundError (gsettings not installed) is suppressed._

- assert detect_system_theme() is None

### TestLinuxTheme.test_linux_gsettings_both_fail_returns_none
_Both gsettings queries failing returns None on Linux._

- assert detect_system_theme() is None

### TestUnknownPlatform.test_unknown_os_returns_none
_Unrecognised OS returns None._

- assert detect_system_theme() is None

### TestDetectSystemDarkMode.test_dark_theme_returns_true
_detect_system_dark_mode returns True when theme is dark._

- assert detect_system_dark_mode() is True

### TestDetectSystemDarkMode.test_light_theme_returns_false
_detect_system_dark_mode returns False when theme is light._

- assert detect_system_dark_mode() is False

### TestDetectSystemDarkMode.test_none_theme_returns_none
_detect_system_dark_mode returns None when theme is unknown._

- assert detect_system_dark_mode() is None

### TestDetectSystemDarkMode.test_winreg_import_error_fallback
_winreg=None is set when the module is unavailable at import time._

- assert su.winreg is None

## tests/unit/test_template_fusing_preview.py

### test_update_template_preview_snaps_points_for_fusing
_Preview geometry snaps nearby polygon vertices to existing scene atoms._

- assert passed_points[0].x() == 0.0
- assert passed_points[0].y() == 0.0
- assert passed_points[1].x() == 10.5
- assert passed_points[1].y() == 0.5

## tests/unit/test_template_preview.py

### test_preview_item_init_defaults
_TemplatePreviewItem initialises with all geometry and flag attributes at their defaults._

- assert item.is_aromatic is False
- assert item.is_user_template is False
- assert item.user_template_points == []
- assert item.user_template_bonds == []
- assert item.user_template_atoms == []

### test_set_geometry_updates_polygon
_set_geometry stores the polygon and clears is_user_template and is_aromatic flags._

- assert not item.polygon.isEmpty()
- assert item.is_aromatic is False
- assert item.is_user_template is False

### test_set_geometry_aromatic_flag
_set_geometry sets is_aromatic=True when the flag is passed._

- assert item.is_aromatic is True

### test_set_user_template_geometry
_set_user_template_geometry stores points, bonds, atoms and sets is_user_template._

- assert item.is_user_template is True
- assert item.user_template_points == pts
- assert item.user_template_bonds == bonds
- assert item.user_template_atoms == atoms
- assert item.is_aromatic is False

### test_bounding_rect_empty_polygon
_boundingRect returns a QRectF when the polygon is empty._

- assert isinstance(rect, QRectF)

### test_bounding_rect_regular_polygon
_boundingRect returns a positive-width rect for a non-empty polygon._

- assert rect.width() > 0

### test_bounding_rect_user_template_with_points
_boundingRect returns positive width and height for a user template with points._

- assert rect.width() > 0
- assert rect.height() > 0

### test_bounding_rect_user_template_no_points
_boundingRect returns a QRectF for a user template with no points._

- assert isinstance(rect, QRectF)

### test_paint_none_painter_returns_early
_paint with a None painter returns without raising._


### test_paint_regular_template_non_aromatic
_paint_regular_template draws a non-aromatic polygon without raising._


### test_paint_regular_template_aromatic
_paint_regular_template draws an aromatic polygon without raising._


### test_paint_regular_template_empty_polygon
_paint_regular_template with an empty polygon does not raise._


### test_paint_user_template_no_points_returns_early
_paint_user_template with no points returns without raising._


### test_paint_user_template_single_bond
_paint_user_template renders a single bond without raising._


### test_paint_user_template_double_bond
_paint_user_template renders a double bond without raising._


### test_paint_user_template_triple_bond
_paint_user_template renders a triple bond without raising._


### test_paint_user_template_non_carbon_atom
_paint_user_template renders non-carbon atom labels without raising._


### test_paint_user_template_bond_info_2_elements
_bond_info with only 2 elements (no order) defaults to single bond._


### test_paint_user_template_out_of_range_indices
_Bond indices beyond point list should be skipped without error._


### test_preview_view_init_defaults
_TemplatePreviewView initialises with original_scene_rect, template_data, and parent_dialog all None._

- assert view.original_scene_rect is None
- assert view.template_data is None
- assert view.parent_dialog is None

### test_preview_view_set_template_data
_set_template_data stores template_data and parent_dialog on the view._

- assert view.template_data == data
- assert view.parent_dialog is parent

### test_preview_view_resize_no_scene_rect_no_timer
_resizeEvent does not schedule a timer when original_scene_rect is None._

- mock_timer.assert_not_called()

### test_preview_view_resize_with_scene_rect_schedules_timer
_resizeEvent schedules a single-shot timer when original_scene_rect is set._

- mock_timer.assert_called_once()

### test_refit_view_no_rect_is_noop
_refit_view does not raise when original_scene_rect is None._


### test_refit_view_empty_rect_is_noop
_refit_view does not raise when original_scene_rect is an empty QRectF._


### test_refit_view_valid_rect
_refit_view does not raise for a valid non-empty QRectF._


### test_show_event_no_rect_no_timer
_showEvent does not schedule a timer when original_scene_rect is None._

- mock_timer.assert_not_called()

### test_show_event_with_rect_schedules_timer
_showEvent schedules a single-shot timer when original_scene_rect is set._

- mock_timer.assert_called_once()

### test_redraw_no_data_returns_early
_redraw_with_current_size returns without raising when template_data is None._


### test_redraw_with_data_calls_draw_template_preview
_redraw_with_current_size calls parent_dialog.draw_template_preview with the stored data._

- parent.draw_template_preview.assert_called_once()

### test_load_template_file_valid
_load_template_file returns the parsed dict for a valid JSON template file._

- assert result == data

### test_load_template_file_missing
_load_template_file returns None for a nonexistent file path._

- assert result is None

### test_load_template_file_invalid_json
_load_template_file returns None when the file contains invalid JSON._

- assert result is None

### test_save_template_file_success
_save_template_file writes JSON to disk and returns True on success._

- assert result is True
- assert json.loads(fp.read_text()) == data

### test_save_template_file_oserror
_save_template_file returns False when the file cannot be written._

- assert result is False

### test_draw_template_preview_no_atoms_adds_placeholder
_draw_template_preview adds a placeholder text item when the template has no atoms._

- assert sc.items()

### test_draw_template_preview_single_bond
_draw_template_preview renders a single bond without raising._

- assert sc.items()

### test_draw_template_preview_double_bond
_draw_template_preview renders a double bond without raising._

- assert sc.items()

### test_draw_template_preview_triple_bond
_draw_template_preview renders a triple bond without raising._

- assert sc.items()

### test_draw_template_preview_non_carbon_atom
_draw_template_preview renders non-carbon atom labels without raising._

- assert sc.items()

### test_draw_template_preview_zero_size_mol
_Single atom (mol_size=0) falls back to scale_factor=1.0 without crash._


### test_convert_structure_to_template
_convert_structure_to_template returns a dict with the correct name, atoms, bonds, and format._

- assert result['name'] == 'TestMol'
- assert len(result['atoms']) == 2
- assert len(result['bonds']) == 1
- assert result['bonds'][0]['order'] == 2
- assert result['format'] == 'PME Template'

### test_cleanup_template_mode_resets_selected
_cleanup_template_mode clears selected_template and disables the delete button._

- assert dlg.selected_template is None
- assert dlg.delete_button.isEnabled() is False

### test_cleanup_template_mode_calls_set_mode
_cleanup_template_mode calls ui_manager.set_mode_and_update_toolbar with 'atom_C'._

- mw.ui_manager.set_mode_and_update_toolbar.assert_called_once_with('atom_C')

### test_cleanup_template_mode_no_ui_manager_no_crash
_cleanup_template_mode does not raise when ui_manager is absent._


### test_cleanup_template_mode_scene_reset
_cleanup_template_mode resets scene.mode and scene.user_template_data._

- assert scene.mode == 'atom_C'
- assert scene.user_template_data is None

### test_select_template_sets_selected_and_enables_delete
_select_template stores the template data and enables the delete button._

- assert dlg.selected_template == data
- assert dlg.delete_button.isEnabled() is True

### test_select_template_activates_template_mode
_select_template calls _activate_template_mode with the chosen template data._

- mock_act.assert_called_once_with(data)

### test_use_template_calls_activate
_use_template calls _activate_template_mode and sets selected_template._

- mock_act.assert_called_once_with(data)
- assert dlg.selected_template == data

### test_fit_preview_view_safely_valid_rect
_fit_preview_view_safely calls fitInView for a non-empty rect._

- view.fitInView.assert_called_once()

### test_fit_preview_view_safely_empty_rect
_fit_preview_view_safely skips fitInView for an empty QRectF._

- view.fitInView.assert_not_called()

### test_fit_preview_view_safely_none_view
_fit_preview_view_safely does not raise when view is None._


### test_refit_all_previews_no_items_no_crash
_refit_all_previews does not raise when the template grid is empty._


### test_refit_all_previews_calls_redraw
_refit_all_previews calls redraw_with_current_size on each TemplatePreviewView child._

- preview_view.redraw_with_current_size.assert_called_once()

### test_delete_selected_no_selection_returns_early
_delete_selected_template returns early without prompting when no template is selected._

- mock_q.assert_not_called()

## tests/unit/test_template_snapping_logic.py

### test_update_template_preview_uses_template_snapping_distance
_update_template_preview calls find_atom_near with template_snapping_distance_2d._

- assert len(scene.find_atom_near_args) == 1
- assert tol == 25.0

### test_update_user_template_preview_uses_template_snapping_distance
_update_user_template_preview calls find_atom_near with template_snapping_distance_2d._

- assert len(scene.find_atom_near_args) == 1
- assert tol == 30.0

### test_keyboard_key_4_uses_template_snapping_distance
_Key 4 (template ring) uses template_snapping_distance_2d for atom snap._

- assert len(scene.find_atom_near_args) == 1
- assert tol == 22.0

### test_keyboard_key_4_uses_template_snapping_when_fusing_disabled
_Key 4 still uses template_snapping_distance_2d even when fusing is disabled._

- assert len(scene.find_atom_near_args) == 1
- assert tol == 22.0

### test_keyboard_other_keys_use_bond_snapping_distance
_Non-template key presses use bond_snapping_distance_2d for atom snap._

- assert len(scene.find_atom_near_args) == 1
- assert tol == 8.0

## tests/unit/test_translation_dialog.py

### TestTabSwitching.test_tab_change_clears_selection
_Switching tabs clears the selected atoms set._

- assert len(dlg.selected_atoms) == 0

### TestTabSwitching.test_tab_change_during_init_is_noop
__on_tab_changed is a no-op while _is_initializing is True._

- assert 0 in dlg.selected_atoms

### TestAbsolutePicking.test_abs_pick_supports_multiple_selection
__abs_on_atom_picked accumulates atoms rather than replacing the selection._

- assert dlg.selected_atoms == {0, 1}

### TestAbsolutePicking.test_abs_pick_populates_inputs
__abs_on_atom_picked fills the absolute coordinate inputs with the atom's position._

- assert dlg.abs_x_input.text() == f'{pos[0][0]:.4f}'
- assert dlg.abs_y_input.text() == f'{pos[0][1]:.4f}'
- assert dlg.abs_z_input.text() == f'{pos[0][2]:.4f}'

### TestAbsoluteHelpers.test_abs_clear_resets_inputs
__abs_clear_selection resets coordinate inputs to 0.000 and empties selected_atoms._

- assert dlg.abs_x_input.text() == '0.000'
- assert dlg.abs_y_input.text() == '0.000'
- assert dlg.abs_z_input.text() == '0.000'
- assert len(dlg.selected_atoms) == 0

### TestAbsoluteHelpers.test_set_origin_fills_zeros
__set_origin fills all absolute coordinate inputs with '0.0000'._

- assert dlg.abs_x_input.text() == '0.0000'
- assert dlg.abs_y_input.text() == '0.0000'
- assert dlg.abs_z_input.text() == '0.0000'

### TestAbsoluteHelpers.test_move_mol_toggled_changes_button_label
__on_move_mol_toggled updates the apply button label based on the checkbox state._

- assert dlg.abs_apply_btn.text() == 'Move Selected'
- assert dlg.abs_apply_btn.text() == 'Move Molecule'

### TestApplyAbsolute.test_no_atom_shows_warning
_apply_absolute shows a warning when no atom is selected._

- mb.warning.assert_called_once()

### TestApplyAbsolute.test_bad_input_shows_warning
_apply_absolute shows a warning when a coordinate input is not a valid number._

- mb.warning.assert_called_once()

### TestApplyAbsolute.test_move_molecule_shifts_all_atoms
_apply_absolute with move_mol=True shifts all atoms by the vector to the target._

- assert after[i] == pytest.approx(before[i] + delta, abs=0.0001)

### TestApplyAbsolute.test_move_selected_only_shifts_selected_atoms
_apply_absolute with move_mol=False shifts only the selected atoms to the target._

- assert after[0] == pytest.approx(before[0] + expected_delta, abs=0.0001)
- assert after[1] == pytest.approx(before[1] + expected_delta, abs=0.0001)
- assert after[i] == pytest.approx(before[i], abs=0.0001)

### TestApplyAbsolute.test_apply_absolute_pushes_undo
_apply_absolute calls push_undo_state on the edit_actions_manager._

- mw.edit_actions_manager.push_undo_state.assert_called()

### TestDeltaPicking.test_delta_pick_adds_atom
__delta_on_atom_picked adds an atom to selected_atoms._

- assert 0 in dlg.selected_atoms

### TestDeltaPicking.test_delta_repick_removes_atom
__delta_on_atom_picked removes an atom that was already selected._

- assert 0 not in dlg.selected_atoms

### TestDeltaPicking.test_delta_picks_accumulate
_Successive _delta_on_atom_picked calls accumulate distinct atoms._

- assert dlg.selected_atoms == {0, 1, 2}

### TestSelectAllAtoms.test_select_all_selects_every_atom
_select_all_atoms selects every atom in the molecule._

- assert dlg.selected_atoms == set(range(mol.GetNumAtoms()))

### TestApplyTranslation.test_no_atoms_shows_warning
_apply_translation shows a warning when selected_atoms is empty._

- mb.warning.assert_called_once()

### TestApplyTranslation.test_bad_input_shows_warning
_apply_translation shows a warning when a delta input is not a valid number._

- mb.warning.assert_called_once()

### TestApplyTranslation.test_zero_vector_is_noop
_apply_translation with a (0,0,0) delta leaves all atom positions unchanged._

- assert after == pytest.approx(before, abs=1e-06)

### TestApplyTranslation.test_delta_shifts_selected_atoms
_apply_translation shifts selected atoms by the given delta, leaving others unchanged._

- assert after[idx] == pytest.approx(before[idx] + [1.0, 2.0, 3.0], abs=0.0001)
- assert after[idx] == pytest.approx(before[idx], abs=0.0001)

### TestApplyTranslation.test_apply_translation_pushes_undo
_apply_translation calls push_undo_state on the edit_actions_manager._

- mw.edit_actions_manager.push_undo_state.assert_called()

### TestUpdateDisplay.test_abs_tab_no_atom_disables_button
_update_display disables the absolute apply button when no atom is selected._

- assert not dlg.abs_apply_btn.isEnabled()

### TestUpdateDisplay.test_abs_tab_one_atom_enables_button
_update_display enables the absolute apply button when one atom is selected._

- assert dlg.abs_apply_btn.isEnabled()

### TestUpdateDisplay.test_delta_tab_no_atom_disables_button
_update_display disables the delta apply button when no atom is selected._

- assert not dlg.apply_button.isEnabled()

### TestUpdateDisplay.test_delta_tab_one_atom_enables_button
_update_display enables the delta apply button when one atom is selected._

- assert dlg.apply_button.isEnabled()

### TestPreselectedAtoms.test_single_preselected_goes_to_abs_tab
_A single preselected atom switches the dialog to the absolute tab._

- assert dlg.tabs.currentIndex() == 0
- assert dlg.selected_atoms == {0}

### TestPreselectedAtoms.test_multiple_preselected_goes_to_delta_tab
_Multiple preselected atoms switch the dialog to the delta tab._

- assert dlg.tabs.currentIndex() == 1
- assert {0, 1, 2}.issubset(dlg.selected_atoms)

## tests/unit/test_ui_manager_robustness.py

### test_close_event_robustness
_Test that handle_close_event handles missing attributes gracefully._

- assert result is True

### test_enable_3d_features_robustness
_Test that _enable_3d_features handles missing widgets._


### test_set_mode_robustness
_Test that set_mode handles template preview visibility._

- assert ui.host.init_manager.scene.mode == 'atom_C'
- assert ui.host.init_manager.scene.template_preview.hide.called

## tests/unit/test_utils_sip.py

### test_sip_isdeleted_safe_valid_obj
_Test safe check with a valid object (mocked)._

- assert sip_isdeleted_safe(obj) is False

### test_sip_isdeleted_safe_deleted_obj
_Test safe check with a deleted object._

- assert sip_isdeleted_safe(obj) is True

### test_sip_isdeleted_safe_exception
_Test safe check when an exception occurs._

- assert sip_isdeleted_safe(obj) is False

### test_sip_isdeleted_safe_no_sip
_Test safe check when _sip_isdeleted is None (sip import failed)._

- assert result is False

## tests/unit/test_view_3d.py

### test_view_3d_draw_standard_3d_style
_Verify that draw_standard_3d_style clears the plotter and constructs the correct VTK meshes._

- view3d.plotter.clear.assert_called()
- view3d.plotter.set_background.assert_called_with('#ffffff')
- assert view3d.plotter.add_mesh.call_count >= 1
- view3d.plotter.render.assert_called()
- assert hasattr(view3d, 'atom_positions_3d')
- assert isinstance(view3d.atom_positions_3d, np.ndarray)
- assert len(view3d.atom_positions_3d) == 3

### test_view_3d_draw_none
_Verify that calling draw with None safely clears the renderer._

- view3d.plotter.clear.assert_called()
- view3d.plotter.render.assert_called()
- assert view3d.current_mol is None

### test_original_id_mode_hides_rdkit_id_labels
_Original-ID labels should not fall back to RDKit index labels._

- view3d.plotter.add_point_labels.assert_called_once()
- assert label_texts == ['42']
- view3d.plotter.add_text.assert_called_once()
- assert view3d.plotter.add_text.call_args.args[0] == 'ID'

### test_fit_to_view_empty
_Verify that fit_to_view resets zoom when the scene has no items._

- view3d.reset_zoom.assert_called_once()

### test_fit_to_view_with_items
_Verify that fit_to_view calculates bounding box and fits view when scene has items._

- view_2d.setTransformationAnchor.assert_called()
- view_2d.setResizeAnchor.assert_called()
- view_2d.fitInView.assert_called_once()

### test_add_3d_atom_glyphs_styles
_Verify radii and resolutions for different styles in _add_3d_atom_glyphs._

- assert mock_pv.PolyData.call_count >= 1
- assert 'radii' in mock_poly.__setitem__.call_args_list[1][0]
- assert np.isclose(rad_array[0], 0.51)
- assert rad_array[0] > 1.5
- assert np.isclose(rad_array[0], 0.15)

### test_add_3d_atom_glyphs_stick_split
_Verify that terminal multiple bonds lead to atom splitting in stick mode._

- assert mock_pv.PolyData.call_count >= 2
- assert len(new_positions) == 4

### test_add_3d_bond_cylinders_basic
_Verify single, double, and triple bond generation._

- assert mock_pv.PolyData.call_count >= 1
- assert len(points) == 14

### test_add_3d_bond_cylinders_styles
_Check style-dependent factors (radius/offset factors) in bond drawing._

- assert np.allclose(radii, 0.08)
- assert np.allclose(radii, 0.09)

### test_add_3d_bond_cylinders_overrides
_Verify plugin bond color overrides._

- assert np.array_equal(colors[0], [255, 0, 0])

### test_add_3d_aromatic_rings
_Verify aromatic torus generation._

- assert mock_pv.Spline.call_count == 1
- assert view3d.plotter.add_mesh.call_count == 1

### test_calculate_double_bond_offset
_Verify neighbor-based plane calculation for double bond offset._

- assert len(offset) == 3
- assert np.isclose(np.linalg.norm(offset), 1.0)

### test_show_ez_labels_3d
_Verify EZ label detection and discrepancy marking._

- assert 'E' in args[1]
- assert '?' in args[1]

### test_chiral_labels_logic
_Verify chiral label toggling and update._

- assert view3d.show_chiral_labels is True
- assert atom_item.chiral_label == 'S'

### test_color_overrides
_Verify color override API functions._

- assert view3d._plugin_bond_color_overrides[0] == '#FF0000'
- view3d.draw_molecule_3d.assert_called()
- assert view3d._plugin_color_overrides[0] == '#00FF00'
- assert view3d.draw_molecule_3d.call_count == 2

## tests/unit/test_worker_robustness.py

### test_worker_halt_logic
_Test that worker correctly identifies halt requests._


### test_iterative_optimize_robustness
_Test that iterative_optimize handles RDKit failures by returning False instead of crashing._

- assert not _iterative_optimize(mol, 'UFF', lambda: False, lambda x: None)

### test_direct_conversion_failure_robustness
_Test that _perform_direct_conversion raises ValueError on empty input._


### test_run_calculation_empty_input
_Test that run_calculation emits error on empty input._

- assert error_handler.called
- assert 'No atoms to convert' in str(args[0][1])

## tests/unit/test_zoomable_view.py

### test_init_panning_state
_ZoomableView initialises with _is_panning=False and zero scroll offsets._

- assert view._is_panning is False
- assert view._pan_start_scroll_h == 0
- assert view._pan_start_scroll_v == 0

### test_init_drag_mode
_ZoomableView initialises with NoDrag drag mode._

- assert view.dragMode() == ZoomableView.DragMode.NoDrag

### test_init_scroll_bars_always_on
_ZoomableView initialises with both scroll bars always visible._

- assert view.verticalScrollBarPolicy() == Qt.ScrollBarPolicy.ScrollBarAlwaysOn
- assert view.horizontalScrollBarPolicy() == Qt.ScrollBarPolicy.ScrollBarAlwaysOn

### test_wheel_ctrl_zoom_in_scales_up
_Ctrl+scroll-up increases the view scale._

- assert view.transform().m11() > before
- ev.accept.assert_called_once()

### test_wheel_ctrl_zoom_out_scales_down
_Ctrl+scroll-down decreases the view scale._

- assert view.transform().m11() < mid

### test_wheel_no_ctrl_passes_to_super
_Scroll without Ctrl passes the event to the parent class and leaves scale unchanged._

- assert view.transform().m11() == pytest.approx(before)

### test_wheel_ctrl_does_not_exceed_max_scale
_Repeated Ctrl+scroll-up clamps the scale near the maximum limit._

- assert view.transform().m11() <= 20.0 * 1.1 + 0.1

### test_wheel_ctrl_does_not_go_below_min_scale
_Repeated Ctrl+scroll-down clamps the scale near the minimum limit._

- assert view.transform().m11() >= 0.05 / 1.1 - 0.01

### test_mouse_press_middle_button_starts_pan
_Middle-button press starts panning and accepts the event._

- assert view._is_panning is True
- ev.accept.assert_called_once()

### test_mouse_press_shift_left_starts_pan
_Shift+left-button press starts panning and accepts the event._

- assert view._is_panning is True
- ev.accept.assert_called_once()

### test_mouse_press_pan_sets_cursor_closed_hand
_Starting a pan changes the cursor to ClosedHandCursor._

- assert view.cursor().shape() == Qt.CursorShape.ClosedHandCursor

### test_mouse_press_left_no_shift_passes_to_super
_Left-button press without Shift does not start panning._

- assert view._is_panning is False

### test_mouse_press_none_event_returns_early
_mousePressEvent with None does not raise and does not start panning._

- assert view._is_panning is False

### test_mouse_move_while_panning_updates_scrollbars
_Mouse move during panning accepts the event and updates scroll positions._

- move_ev.accept.assert_called_once()

### test_mouse_move_not_panning_passes_to_super
_Mouse move without active pan passes the event to the parent class._

- assert view._is_panning is False

### test_mouse_move_none_event_returns_early
_mouseMoveEvent with None does not raise._


### test_mouse_release_ends_pan
_Mouse release while panning ends the pan and accepts the event._

- assert view._is_panning is True
- assert view._is_panning is False
- rel_ev.accept.assert_called_once()

### test_mouse_release_restores_arrow_cursor_in_select_mode
_Ending a pan in select mode restores the ArrowCursor._

- assert view.cursor().shape() == Qt.CursorShape.ArrowCursor

### test_mouse_release_restores_cross_cursor_in_atom_mode
_Ending a pan in atom mode restores the CrossCursor._

- assert view.cursor().shape() == Qt.CursorShape.CrossCursor

### test_mouse_release_restores_cross_cursor_in_bond_mode
_Ending a pan in bond mode restores the CrossCursor._

- assert view.cursor().shape() == Qt.CursorShape.CrossCursor

### test_mouse_release_restores_cross_cursor_in_charge_mode
_Ending a pan in charge mode restores the CrossCursor._

- assert view.cursor().shape() == Qt.CursorShape.CrossCursor

### test_mouse_release_restores_arrow_for_unknown_mode
_Ending a pan in an unknown scene mode restores the ArrowCursor._

- assert view.cursor().shape() == Qt.CursorShape.ArrowCursor

### test_mouse_release_not_panning_passes_to_super
_Mouse release when not panning passes the event to the parent class._

- assert view._is_panning is False
- assert view._is_panning is False

### test_viewport_event_non_native_gesture_passes_to_super
_viewportEvent with a non-NativeGesture event passes to the parent class._

- assert isinstance(result, bool)

### test_viewport_event_zoom_native_gesture_scales
_viewportEvent with ZoomNativeGesture scales the view and returns True._

- assert result is True
- assert view.transform().m11() > before

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

### test_calculation_worker_opt_failure_emits_error
_Test that optimization failure emits an error string instead of failing silently or hanging._

- assert mock_mmff_props.called
- assert res_mol is not None

### test_calculation_worker_mmff_variants
_Test switching between MMFF94 and MMFF94s._

- assert found

### test_calculation_worker_obabel_fallback_mocked
_Test the Open Babel fallback path by mocking availability and pybel._

- assert mock_run.called
- assert sys.executable in args[0]

### test_calculation_worker_complex_direct_h_placement
_Test direct mode with 4 hydrogens on a carbon to hit the rotation/offset logic._

- assert mol_3d.GetNumAtoms() >= 5

### test_worker_halt_error_is_exception
_WorkerHaltError is a proper Exception subclass._

- assert isinstance(err, Exception)
- assert str(err) == 'test halt'

### test_worker_halt_error_not_caught_by_generic
_WorkerHaltError should propagate through except Exception if re-raised._


### test_iterative_optimize_mmff
__iterative_optimize: MMFF method converges on ethane._

- assert result is True

### test_iterative_optimize_uff
__iterative_optimize: UFF method converges on ethane._

- assert result is True

### test_iterative_optimize_mmff94_variant
__iterative_optimize: MMFF94 (non-s) variant._

- assert result is True

### test_iterative_optimize_unknown_method
__iterative_optimize: unknown method returns False._

- assert result is False

### test_iterative_optimize_halt_during_optimization
__iterative_optimize: raises WorkerHaltError when halted during chunks._


### test_iterative_optimize_props_none
__iterative_optimize: returns False when MMFF properties cannot be computed._

- assert result is False

### test_iterative_optimize_ff_none
__iterative_optimize: returns False when UFF force field is None._

- assert result is False

### test_adjust_collision_avoidance_single_fragment
__adjust_collision_avoidance: single fragment returns immediately (no-op)._

- assert len(status_msgs) == 0

### test_adjust_collision_avoidance_multi_fragment
__adjust_collision_avoidance: two overlapping fragments get separated._

- assert len(frags) == 2
- assert len(status_msgs) >= 2
- assert 'Resolving' in status_msgs[0]
- assert 'completed' in status_msgs[-1]

### test_adjust_collision_avoidance_halt
__adjust_collision_avoidance: raises WorkerHaltError when halted._


### test_calculation_worker_direct_with_optimize
_Integration test: Direct conversion mode with do_optimize=True._

- assert mol_3d is not None
- assert mol_3d.GetNumConformers() >= 1

### test_calculation_worker_optimized_result_better
_Integration test: Optimized ethane should have a C-C bond length_

- assert 1.3 < dist < 1.7

### test_calculation_worker_optimize_only_mmff
_Integration test: optimize_only mode with MMFF method._

- assert mol_3d is not None
- assert mol_3d.GetNumConformers() >= 1

### test_calculation_worker_optimize_only_uff
_Integration test: optimize_only mode with UFF method._

- assert mol_3d is not None

### test_calculation_worker_optimize_only_default
_Integration test: optimize_only mode with default (no optimization_method)._

- assert mol_3d is not None

### test_calculation_worker_optimize_only_mmff94_variant
_Integration test: optimize_only mode with MMFF94 (non-s) variant._

- assert mol_3d is not None

### test_calculation_worker_status_signals
_Integration test: verify status_update signals are emitted during conversion._

- assert len(status_messages) >= 1
- assert any(('3D' in msg or 'Creating' in msg for msg in status_messages))

### test_calculation_worker_multi_fragment_rdkit
_Integration test: multi-fragment molecule triggers collision avoidance in RDKit path._

- assert mol_3d is not None
- assert len(collision_msgs) >= 1

### test_calculation_worker_direct_dash_stereo
_Integration test: direct mode with a dash bond._

- assert p3.z == pytest.approx(-1.5)

### test_calculation_worker_direct_mmff94_rdkit_variant
_Integration test: direct mode with do_optimize=True and MMFF94_RDKIT variant._

- assert mol_3d is not None
- assert mol_3d.GetNumConformers() >= 1

### test_calculation_worker_fallback_to_direct_no_obabel
_Integration test: fallback mode with OBABEL_AVAILABLE=False and_

- assert mol_3d is not None
- assert mol_3d.GetNumConformers() >= 1
- assert any(('direct' in m.lower() or 'Direct' in m for m in status_messages))

### test_calculation_worker_fallback_to_direct_with_optimize
_Integration test: fallback-to-direct with do_optimize=True._

- assert mol_3d is not None
- assert mol_3d.GetNumConformers() >= 1

## tests/integration/test_chiral_labels_integration.py

### test_chiral_labels_toggle_3d
_Test that toggling 'Show Chiral Labels' correctly displays/hides labels in 3D._

- assert window.view_3d_manager.show_chiral_labels is False
- assert not labels_drawn()
- assert window.view_3d_manager.show_chiral_labels is True
- assert len(chiral_centers) == 1
- assert initial_label in ['R', 'S']
- assert len(chiral_call.args) > 1
- assert initial_label in labels

### test_chiral_labels_mirror_inversion_3d
_Test that mirror transformation inverts the chiral label in 3D._

- assert initial_label in ['R', 'S']
- assert len(chiral_centers) == 1
- assert new_label != initial_label
- assert new_label in ['R', 'S']
- assert len(chiral_call.args) > 1
- assert new_label in labels
- assert initial_label not in labels

## tests/integration/test_headless_install.py

### test_headless_install_success
_Test successful headless plugin installation._

- assert 'PLUGIN INSTALLATION (HEADLESS)' in captured.out
- assert 'CLI Test Plugin' in captured.out
- assert 'Success:' in captured.out
- assert cm.value.code == 0

### test_headless_install_abort
_Test aborted headless plugin installation (answering 'n')._

- mock_pm.install_plugin.assert_not_called()
- assert 'Installation aborted.' in captured.out
- assert cm.value.code == 0

### test_headless_install_invalid_path
_Test error message when plugin path is invalid._

- assert 'Error: Plugin path not found' in captured.out
- assert cm.value.code == 1

## tests/integration/test_plugin_menu_rebuild.py

### TestRebuildCompleteness.test_all_six_steps_run_when_all_registries_empty
_Even with empty registries, all six rebuild methods execute._

- assert called == ['menu_actions', 'toolbar_actions', 'export_actions', 'file_openers', 'analysis_tools', 'style_menu']

### TestRebuildCompleteness.test_all_six_steps_run_when_all_registries_populated
_With all registries populated, all six steps still run._

- assert set(called) == {'menu_actions', 'toolbar_actions', 'export_actions', 'file_openers', 'analysis_tools', 'style_menu'}

### TestRebuildCompleteness.test_separator_flag_reset_before_rebuild
_rebuild_plugin_menus resets the separator flag to False before running._

- assert im.plugin_menubar_separator_added is False

### TestRebuildCompleteness.test_step_exception_does_not_abort_remaining_steps
_If one step raises, all subsequent steps still execute._

- assert 'toolbar' in survived
- assert 'export' in survived
- assert 'style' in survived

### TestPluginMenuManagerRouting.test_main_window_proxy_returns_init_manager_pmm
_MainWindow.plugin_menu_manager property reads from init_manager._

- assert prop is not None and isinstance(prop, property)
- assert result is pmm_mock

### TestPluginMenuManagerRouting.test_plugin_manager_rebuild_calls_main_window_pmm
_PluginManager.rebuild_plugin_menus() calls main_window.plugin_menu_manager._

- pmm_mock.rebuild_plugin_menus.assert_called_once()

### TestPluginMenuManagerRouting.test_plugin_manager_rebuild_no_main_window_does_nothing
_PluginManager.rebuild_plugin_menus() with no main_window is a no-op._


### TestPluginMenuManagerRouting.test_main_init_manager_has_no_wrapper_methods
_MainInitManager no longer exposes the old wrapper methods directly._

- assert not hasattr(MainInitManager, method)

### TestUpdatePluginMenuIntegration.test_discover_plugins_called_with_host
_update_plugin_menu calls discover_plugins with the init_manager host._

- pm.discover_plugins.assert_called_once_with(im.host)

### TestUpdatePluginMenuIntegration.test_calls_seven_integration_methods
_update_plugin_menu invokes all seven integration methods exactly once._


### TestUpdatePluginMenuIntegration.test_does_nothing_when_plugin_manager_is_none
_update_plugin_menu is a no-op when plugin_manager is None._

- plugin_menu.clear.assert_not_called()

### TestExportActionEndToEnd.test_export_action_appears_in_button_menu
_A registered export action is added to the export button menu._

- export_menu.addSeparator.assert_called_once()
- export_menu.addAction.assert_called_once()
- assert added_action.text() == 'Export XYZ'

### TestExportActionEndToEnd.test_multiple_export_actions_all_added
_All registered export actions are each added to the export menu._

- assert im.export_button.menu.return_value.addAction.call_count == 2

### TestAnalysisToolEndToEnd.test_analysis_tools_appear_in_analysis_menu
_Registered analysis tools appear as actions in the Analysis menu after rebuild._

- analysis_menu.addSeparator.assert_called_once()
- assert analysis_menu.addAction.call_count == 1
- assert 'NMR Predict' in added.text()
- assert 'NMRPro' in added.text()

## tests/integration/test_trigger_conversion_plugin_wrap.py

### TestTriggerConversionTempMode.test_temp_mode_stored_and_consumed
__trigger_conversion_with_temp_mode sets _pending_conversion_mode_

- assert mock_timer.called
- assert compute.__dict__.get('_pending_conversion_mode') == 'rdkit'

### TestTriggerConversionTempMode.test_pending_mode_consumed_after_trigger_conversion
_After trigger_conversion runs, _pending_conversion_mode is None._

- assert compute.__dict__.get('_pending_conversion_mode') is None

### TestTriggerConversionTempMode.test_pending_mode_takes_priority_over_settings
__pending_conversion_mode overrides 3d_conversion_mode setting._

- assert options['conversion_mode'] == 'fallback'

### TestTriggerConversionTempMode.test_settings_mode_used_when_no_pending_mode
_When _pending_conversion_mode is absent, the settings value is used._

- assert options['conversion_mode'] == 'obabel'

### TestPluginWrappedTriggerConversion.test_no_typeerror_when_plugin_wraps_without_kwargs
_Calling _trigger_conversion_with_temp_mode must not raise TypeError_


### TestPluginWrappedTriggerConversion.test_mode_reaches_worker_through_plugin_wrapper
_The conversion mode set via _trigger_conversion_with_temp_mode_

- assert hasattr(compute, '_captured_options')
- assert compute._captured_options['conversion_mode'] == 'fallback'

## tests/gui/test_additional_dialogs_launch.py

### test_planarize_dialog_launch
_PlanarizeDialog opens with the correct window title._

- assert dialog.windowTitle() == 'Planarize'

### test_mirror_dialog_launch
_MirrorDialog opens with the correct window title._

- assert dialog.windowTitle() == 'Mirror Molecule'

### test_align_plane_dialog_launch
_AlignPlaneDialog opens with a plane-specific window title._

- assert dialog.windowTitle() == f'Align to {plane_names[plane]} Plane'

### test_constrained_optimization_dialog_launch
_ConstrainedOptimizationDialog opens with the correct window title._

- assert dialog.windowTitle() == 'Constrained Optimization'

### test_periodic_table_dialog_launch
_PeriodicTableDialog opens with the correct window title._

- assert dialog.windowTitle() == 'Select an Element'

## tests/gui/test_dialog_launch.py

### test_bond_length_dialog_launch
_BondLengthDialog opens with the correct window title._

- assert dialog.windowTitle() == 'Adjust Bond Length'

### test_angle_dialog_launch
_AngleDialog opens with the correct window title._

- assert dialog.windowTitle() == 'Adjust Angle'

### test_dihedral_dialog_launch
_DihedralDialog opens with the correct window title._

- assert dialog.windowTitle() == 'Adjust Dihedral Angle'

### test_alignment_dialog_launch
_AlignmentDialog opens with axis-specific window title._

- assert dialog.windowTitle() == 'Align to X-axis'

### test_translation_dialog_launch
_TranslationDialog opens with the correct window title._

- assert dialog.windowTitle() == 'Translate Atoms'

### test_move_group_dialog_launch
_MoveGroupDialog opens with the correct window title._

- assert dialog.windowTitle() == 'Move Group'

## tests/gui/test_edit_actions_gui_extended.py

### test_rotate_2d_dialog_init
_Test Rotate2DDialog GUI initialization in the GUI test environment._

- assert dialog.windowTitle() == 'Rotate 2D'
- assert dialog.angle_spin.value() == 45
- assert dialog.slider.value() == 45
- assert dialog.slider.value() == 90
- assert dialog.angle_spin.value() == -30
- assert dialog.get_angle() == -30.0

## tests/gui/test_main_app.py

### test_app_launch
_MainWindow: Verify application launches correctly._

- assert window.isVisible()
- assert 'MoleditPy Ver.' in window.windowTitle()

### test_mode_change_atom
_Toolbar: Verify mode changes upon clicking atom buttons._

- assert scene.mode == 'atom_C'
- assert scene.mode == 'atom_N'
- assert scene.current_atom_symbol == 'N'
- assert window.statusBar().currentMessage() == 'Mode: Draw Atom (N)'

### test_mode_change_bond
_Toolbar: Verify mode changes upon clicking bond buttons._

- assert scene.mode == 'bond_2_0'
- assert scene.bond_order == 2
- assert scene.bond_stereo == 0
- assert window.statusBar().currentMessage() == 'Mode: Draw Bond (Order: 2)'

### test_draw_atom_on_click
_MoleculeScene: Test for drawing an atom upon clicking._

- assert len(window.state_manager.data.atoms) == 1
- assert window.state_manager.data.atoms[atom_id]['symbol'] == 'N'
- assert window.init_manager.scene.atom_items[atom_id].pos() == click_pos

### test_draw_bond_on_drag
_MoleculeScene: Test for drawing a bond upon dragging._

- assert len(window.state_manager.data.atoms) == 2
- assert len(window.state_manager.data.bonds) == 1
- assert (id1, id2) in window.state_manager.data.bonds
- assert window.state_manager.data.bonds[id1, id2]['order'] == 1

### test_draw_bond_to_existing_atom
_MoleculeScene: Test for dragging a bond to an existing atom._

- assert len(window.state_manager.data.atoms) == 2
- assert len(window.state_manager.data.bonds) == 0
- assert len(window.state_manager.data.atoms) == 2
- assert len(window.state_manager.data.bonds) == 1
- assert (0, 1) in window.state_manager.data.bonds

### test_change_atom_symbol_on_click
_MoleculeScene: Test for changing the element symbol of an existing atom._

- assert window.state_manager.data.atoms[0]['symbol'] == 'C'
- assert window.state_manager.data.atoms[0]['symbol'] == 'O'
- assert len(window.state_manager.data.atoms) == 1

### test_change_bond_order_on_click
_MoleculeScene: Test for changing the order of an existing bond._

- assert window.state_manager.data.bonds[0, 1]['order'] == 1
- assert window.state_manager.data.bonds[0, 1]['order'] == 2

### test_delete_atom_on_right_click
_MoleculeScene: Test for deleting an atom upon right-clicking._

- assert len(window.state_manager.data.atoms) == 1
- assert len(window.state_manager.data.atoms) == 0

### test_charge_mode_click
_MoleculeScene: Test for clicking in charge mode._

- assert window.state_manager.data.atoms[0]['charge'] == 0
- assert window.state_manager.data.atoms[0]['charge'] == 1
- assert window.state_manager.data.atoms[0]['charge'] == -1

### test_2d_to_3d_conversion
_2D->3D Conversion: Test for the conversion button._

- assert len(window.state_manager.data.atoms) == 2
- assert len(window.state_manager.data.bonds) == 1
- assert convert_button.isEnabled()
- assert window.view_3d_manager.current_mol is not None
- assert window.init_manager.optimize_3d_button.isEnabled()
- assert window.init_manager.export_button.isEnabled()
- assert window.init_manager.analysis_action.isEnabled()

### test_optimize_3d
_3D Optimization: Test for the 3D optimization button._

- assert window.view_3d_manager.current_mol is not None
- assert window.init_manager.optimize_3d_button.isEnabled()
- assert 'Optimization completed' in msg or 'optimization successful' in msg or 'Process completed' in msg

### test_change_3d_style
_3D Style Change: Test for the style menu._

- assert window.view_3d_manager.current_3d_style == 'ball_and_stick'
- assert style_button is not None
- assert cpk_action is not None
- assert window.view_3d_manager.current_3d_style == 'cpk'

### test_undo_redo
_Undo/Redo: Test for editing operations._

- assert len(window.state_manager.data.atoms) == 0
- assert len(window.state_manager.data.atoms) == 1
- assert len(window.edit_actions_manager.undo_stack) == initial_stack_size + 1
- assert window.init_manager.undo_action.isEnabled() is True
- assert window.init_manager.redo_action.isEnabled() is False
- assert len(window.state_manager.data.atoms) == 0
- assert window.init_manager.redo_action.isEnabled() is True
- assert len(window.state_manager.data.atoms) == 1
- assert window.init_manager.undo_action.isEnabled() is True
- assert window.init_manager.redo_action.isEnabled() is False

### test_clear_all
_Clear All: Test for clearing the entire scene._

- assert len(window.state_manager.data.atoms) == 0
- assert len(window.state_manager.data.bonds) == 0
- assert window.view_3d_manager.current_mol is None
- assert window.host.state_manager.has_unsaved_changes is False
- assert len(window.edit_actions_manager.undo_stack) == 1

### test_copy_paste
_Edit: Test for copy & paste operations._

- assert len(window.state_manager.data.atoms) == 2
- assert len(window.state_manager.data.bonds) == 1
- assert len(scene.selectedItems()) == 3
- assert mime_data.hasFormat(moleditpy.CLIPBOARD_MIME_TYPE)
- assert len(window.state_manager.data.atoms) == 4
- assert len(window.state_manager.data.bonds) == 2
- assert window.state_manager.data.atoms[2]['pos'][0] > 50
- assert window.state_manager.data.atoms[2]['pos'][1] > 0
- assert window.state_manager.data.atoms[3]['pos'][0] > 50
- assert window.state_manager.data.atoms[3]['pos'][1] > 0

### test_file_import_smiles
_File: Test for SMILES import._

- assert len(window.state_manager.data.atoms) == 3
- assert len(window.state_manager.data.bonds) == 2
- assert symbols.count('C') == 2
- assert symbols.count('O') == 1
- assert 'Successfully loaded from SMILES' in window.statusBar().currentMessage()

### test_key_press_change_atom
_Keyboard Shortcut: Change atom symbol via 'O' key._

- assert window.state_manager.data.atoms[0]['symbol'] == 'C'
- assert window.state_manager.data.atoms[0]['symbol'] == 'O'

### test_key_press_change_bond
_Keyboard Shortcut: Change bond order via '2' key._

- assert window.state_manager.data.bonds[0, 1]['order'] == 1
- assert window.state_manager.data.bonds[0, 1]['order'] == 2

### test_radical_mode_toggle
_MoleculeScene: Click test in radical mode._

- assert window.state_manager.data.atoms[0]['radical'] == 0
- assert window.state_manager.data.atoms[0]['radical'] == 1
- assert window.state_manager.data.atoms[0]['radical'] == 2
- assert window.state_manager.data.atoms[0]['radical'] == 0

### test_delete_key_selection
_MoleculeScene: Delete selected items via Delete key._

- assert len(window.state_manager.data.atoms) == 1
- assert len(scene.selectedItems()) == 1
- assert len(window.state_manager.data.atoms) == 0

### test_draw_benzene_template
_MoleculeScene: Draw benzene template._

- assert scene.mode == 'template_benzene'
- assert len(window.state_manager.data.atoms) == 6
- assert len(window.state_manager.data.bonds) == 6
- assert orders == [1, 1, 1, 2, 2, 2]

### test_open_settings_dialog
_MainWindow: Test for opening the settings dialog._

- QDialog.exec.assert_called()

### test_toggle_measurement_mode
_MainWindow: Test for toggling 3D measurement mode._

- assert window.measurement_mode is False
- assert measurement_action is not None
- assert window.measurement_mode is True
- assert window.statusBar().currentMessage().startswith('Measurement mode enabled')
- assert window.measurement_mode is False
- assert window.statusBar().currentMessage() == 'Measurement mode disabled.'

### test_toggle_3d_edit_mode
_MainWindow: Test for toggling 3D drag mode._

- assert window.is_3d_edit_mode is False
- assert edit_3d_action is not None
- assert window.is_3d_edit_mode is True
- assert window.statusBar().currentMessage() == '3D Drag Mode: ON.'
- assert window.is_3d_edit_mode is False
- assert window.statusBar().currentMessage() == '3D Drag Mode: OFF.'

### test_add_remove_hydrogens
_Edit: Add/Remove hydrogens menu test._

- assert len(window.state_manager.data.atoms) == 1
- assert len(window.state_manager.data.atoms) == 5
- assert len(window.state_manager.data.bonds) == 4
- assert symbols.count('H') == 4
- assert len(window.state_manager.data.atoms) == 1
- assert len(window.state_manager.data.bonds) == 0
- assert window.state_manager.data.atoms[0]['symbol'] == 'C'

### test_2d_cleanup
_2D Cleanup: Verify coordinates change upon button click._

- assert not (pos0_before == pos0_after and pos1_before == pos1_after)
- assert '2D structure optimization successful' in window.statusBar().currentMessage()

### test_3d_viewer_mode_mol
_3D Viewer Mode: Integration test for MOL file loading and UI state transition._

- assert window.view_3d_manager.current_mol is not None
- assert window.view_3d_manager.current_mol.GetNumAtoms() == 1
- assert window.ui_manager.is_2d_editable is False
- assert window.init_manager.cleanup_button.isEnabled() is False
- assert get_button(window.init_manager.toolbar, 'N (n)').isEnabled() is False
- assert window.init_manager.optimize_3d_button.isEnabled() is True
- assert window.init_manager.export_button.isEnabled() is True
- assert window.init_manager.analysis_action.isEnabled() is True

### test_open_3d_edit_dialogs
_3D Edit: Verify 3D editing dialogs launch correctly._

- assert window.view_3d_manager.current_mol is not None
- assert window.translation_action.isEnabled() is True
- assert window.align_menu.isEnabled() is True
- assert window.planarize_action.isEnabled() is True
- QDialog.show.assert_called()
- assert len(window.active_3d_dialogs) == 0
- QDialog.show.assert_called()

### test_save_project_as
_Project Save: Test for "Save Project As..."._

- mocker_json_dump.assert_called_once()
- assert window.state_manager.has_unsaved_changes is False
- assert window.init_manager.current_file_path == '/fake/save.pmeprj'
- assert 'Project saved to' in window.statusBar().currentMessage()

### test_open_project
_Project Load: Test for "Open Project..." (.pmeprj)._

- assert len(window.state_manager.data.atoms) == 2
- assert len(window.state_manager.data.bonds) == 1
- assert 0 in window.state_manager.data.atoms
- assert 1 in window.state_manager.data.atoms
- assert window.state_manager.data.atoms[0]['symbol'] == 'C'
- assert (0, 1) in window.state_manager.data.bonds
- assert window.view_3d_manager.current_mol is None
- assert 'Project loaded from' in window.statusBar().currentMessage()

### test_toggle_3d_atom_info
_3D Atom Info Display: Test for toggling ID, Coordinates, and Symbol display._

- assert window.view_3d_manager.current_mol is not None
- assert window.atom_info_display_mode == 'original_id'
- mock_add_labels.assert_called()
- assert window.current_atom_info_labels is not None
- assert window.atom_info_display_mode == 'coords'
- mock_add_labels.assert_called()
- assert window.atom_info_display_mode is None
- assert window.current_atom_info_labels is None

### test_user_template_dialog_save_and_use
_User Templates: Test for opening dialog, saving current structure, and using it._

- assert action_open_dialog is not None
- QDialog.show.assert_called()
- assert len(window.state_manager.data.atoms) == 1
- mocker_json_dump.assert_called_once()
- assert QMessageBox.information.called
- assert any(("Template 'test' saved" in str(a) for a in called_args))
- assert len(window.state_manager.data.atoms) == 2
- assert 1 in window.state_manager.data.atoms
- assert window.state_manager.data.atoms[1]['symbol'] == 'C'
- assert window.state_manager.data.atoms[1]['pos'][0] > 50

### test_implicit_hydrogens_update
_Implicit Hydrogens: Test for automatic updates after drawing operations._

- assert len(window.state_manager.data.atoms) == 1
- assert atom_item.implicit_h_count == 4
- assert len(window.state_manager.data.atoms) == 2
- assert len(window.state_manager.data.bonds) == 1
- assert window.init_manager.scene.atom_items[0].implicit_h_count == 3
- assert window.init_manager.scene.atom_items[1].implicit_h_count == 3

### test_drag_drop_mol_file_on_3d_view
_D&D: Drag & Drop of .mol file onto 3D view area (Mock)._

- mock_load_3d.assert_called_once_with(file_path='/fake/drop.mol')
- mock_event.acceptProposedAction.assert_called_once()

### test_drag_drop_mol_file_on_2d_view
_D&D: Drag & Drop of .mol file onto 2D view area (Mock)._

- mock_load_2d.assert_called_once_with(file_path='/fake/drop.mol')
- mock_event.acceptProposedAction.assert_called_once()

### test_project_save_load_round_trip
_Project Save/Load Round-trip: Integration test to verify structure recovery after saving and reloading._

- assert len(window.state_manager.data.atoms) == 2
- assert window.state_manager.data.atoms[id1]['symbol'] == 'N'
- assert (0, 1) in window.state_manager.data.bonds
- assert save_file.exists()
- assert window.host.state_manager.has_unsaved_changes is False
- assert len(window.state_manager.data.atoms) == 0
- assert len(window.state_manager.data.atoms) == 2
- assert len(window.state_manager.data.bonds) == 1
- assert symbols == ['C', 'N']
- assert bond_data['order'] == 1
- assert bond_data['stereo'] == 0

### test_file_import_smiles_error
_SMILES Import: Test for error handling when invalid SMILES is entered._

- assert status_msg.startswith('Invalid SMILES:')

### test_undo_redo_boundary
_Undo/Redo: Test for stack boundary (operations on empty stack)._

- assert len(window.edit_actions_manager.undo_stack) == 1
- assert len(window.state_manager.data.atoms) == 0

### test_import_invalid_mol_file
_File Load: Error handling for corrupted MOL files._

- assert 'Invalid MOL file format:' in status_msg or 'Error loading file:' in status_msg
- assert len(window.state_manager.data.atoms) == 0

### test_clear_2d_editor_cancel
_Clear All: Test for cancellation in confirmation dialog._

- assert len(window.state_manager.data.atoms) == 1
- assert window.host.state_manager.has_unsaved_changes is True

### test_clipboard_copy_empty_selection
_Copy: Safety test for copy operation with empty selection._

- assert cb.text() == 'initial_text'

## tests/gui/test_main_window_settings.py

### test_load_command_line_file_with_plugin
_Test that load_command_line_file uses plugin openers when available._

- mock_callback.assert_called_once_with(test_file)
- assert window.init_manager.current_file_path == test_file

### test_load_command_line_file_default_extensions
_Test that load_command_line_file handles standard extensions._

- window.io_manager.load_mol_file_for_3d_viewing.assert_called_once_with('test.mol')
- window.io_manager.load_xyz_for_3d_viewing.assert_called_once_with('test.xyz')
- window.io_manager.open_project_file.assert_called_once_with(file_path='test.pmeprj')

### test_update_cpk_colors_from_settings
_Test that CPK colors are updated correctly from settings overrides._

- assert constants.CPK_COLORS['C'] == QColor('#FF0000')
- assert constants.CPK_COLORS_PV['C'] == [1.0, 0.0, 0.0]
- assert constants.CPK_COLORS['C'] == DEFAULT_CPK_COLORS.get('C', constants.CPK_COLORS['C'])

### test_apply_initial_settings
_Test that apply_initial_settings updates scene background and style._

- assert window.init_manager.scene.backgroundBrush().color() == expected_color
- window.view_3d_manager.plotter.set_background.assert_called_with('#112233')

### test_reset_all_settings_flow
_Test the complete settings reset flow._

- window.init_manager._perform_settings_reset.assert_called_once()
- window.init_manager._refresh_ui_after_reset.assert_called_once()

### test_perform_settings_reset_logic
_Test the low-level settings reset (file deletion and reload)._

- assert not settings_file.exists()
- window.init_manager.load_settings.assert_called_once()
- assert window.init_manager.settings_dirty is True

## tests/gui/test_main_window_ui_integration.py

### test_plugin_menu_actions_population
_Test that menu actions registered by plugins are correctly added to the menu bar._

- assert plugins_action is not None
- assert plugins_menu is not None
- assert sub_menu_action is not None
- assert sub_menu is not None
- assert test_action is not None
- assert test_action.shortcut().toString() == 'Ctrl+Shift+P'
- assert plugin_analysis_action is not None

### test_plugin_toolbar_actions_visibility
_Test that the plugin toolbar is shown/hidden and populated correctly._

- assert window.init_manager.plugin_toolbar.isHidden()
- assert not window.init_manager.plugin_toolbar.isHidden()
- assert len(toolbar_actions) == 1
- assert toolbar_actions[0].text() == 'ToolBtn'
- assert toolbar_actions[0].toolTip() == 'Hint'

### test_ui_sync_after_reset
_Test that UI elements (background, checked states) sync correctly after settings reset._

- assert window.init_manager.scene is not None
- assert actual_bg == '#0000FF'
- assert window.opt3d_actions['MMFF_RDKIT'].isChecked()
- assert window.conv_actions['rdkit'].isChecked()

### test_custom_3d_style_integration
_Test that custom 3D styles from plugins appear in the style menu._

- assert style_action is not None
- assert style_action.actionGroup() is not None

### test_integrate_plugin_export_actions
_Test that export actions are added to File/Export menu._

- assert file_action is not None
- assert export_action is not None
- assert any(('Export to Fax' in a.text() for a in btn_menu.actions()))
- assert any(('Export to Fax' in a.text() for a in export_menu.actions()))

### test_integrate_plugin_analysis_tools
_Test integration into the Analysis menu._

- assert analysis_action is not None
- assert found

### test_integrate_plugin_file_openers_ui
_Test integration of plugin openers into the Import menu._

- assert import_action is not None
- assert 'Import .fake' in import_action.text()
- cb.assert_called_once_with('data.fake')
- assert window.init_manager.current_file_path == 'data.fake'

## tests/gui/test_molecule_scene_events.py

### test_bond_stereo_toggle_keys
_Test Z and E keys toggle double bond stereochemistry._

- assert data.bonds[bond_key]['stereo'] == 3
- assert bond_item.stereo == 3
- assert data.bonds[bond_key]['stereo'] == 4
- assert bond_item.stereo == 4

### test_atom_addition_keys
_Test 1, 2, 3 keys add atoms/bonds from selected atom._

- assert len(data.atoms) == 2
- assert len(data.bonds) == 1
- assert data.bonds[bond_key]['order'] == 1
- assert len(other_atom_ids) == 1
- assert len(data.atoms) == 3
- assert len(data.bonds) == 2
- assert b_data['order'] == 2

### test_delete_items_keys
_Test Delete and Backspace keys remove selected items._

- assert a1_id not in data.atoms
- assert a2_id in data.atoms
- assert a2_id not in data.atoms
- assert len(data.atoms) == 0

### test_temp_line_cancellation
_Test Delete key cancels an active temp_line (bond drawing)._

- assert scene.temp_line is None

### test_bonding_to_existing_atom
_Test that pressing 1, 2, 3 bonds to an existing atom if it's nearby._

- assert len(data.atoms) == 2
- assert len(data.bonds) == 1
- assert bond_key in data.bonds
