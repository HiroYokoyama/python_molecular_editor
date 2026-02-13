# Test Assertions Catalog

## tests/unit/test_molecular_data.py

### test_add_atom_returns_incrementing_ids
_Verify atom IDs are sequential integers starting from 0._

- assert (id0, id1, id2) == (0, 1, 2)
- assert len(data.atoms) == 3

### test_add_atom_stores_properties
_Verify all atom properties (symbol, position, charge, radical) are stored._

- assert atom['symbol'] == "N"
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

- assert key == (a, b)  # Should be normalized
- assert status == 'created'

### test_add_bond_preserves_stereo_key_direction
_Stereo bonds should preserve the original key direction._

- assert key == (b, a)  # NOT normalized — direction matters

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

- assert pos.x == pytest.approx(px * ANGSTROM_PER_PIXEL, abs=1e-6)
- assert pos.y == pytest.approx(-py * ANGSTROM_PER_PIXEL, abs=1e-6)  # Y is inverted

### test_to_rdkit_mol_preserves_original_atom_id
_Verify _original_atom_id property is set on RDKit atoms._

- assert atom.HasProp("_original_atom_id")
- assert orig_id in [id0, id1]

### test_to_rdkit_mol_empty_returns_none
_Empty MolecularData should return None._

- assert data.to_rdkit_mol() is None

### test_to_mol_block_contains_all_atoms
_Verify MOL block text contains all expected atoms._

- assert mol_block is not None
- assert "V2000" in mol_block or "MoleditPy" in mol_block
- assert symbols == {"C", "N", "O"}

### test_molecular_weight_matches_rdkit
_Build acetic acid and compare MW against RDKit reference._

- assert Descriptors.HeavyAtomMolWt(mol) == pytest.approx(Descriptors.HeavyAtomMolWt(ref), abs=0.01)

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
- assert all(a['symbol'] == 'C' for a in data.atoms.values())

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

- assert "Invalid" in last_msg or "Error" in last_msg

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

- assert "Invalid" in last_msg or "Error" in last_msg

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

## tests/unit/test_plugin_manager.py

### test_plugin_info_extracts_all_fields
_Verify AST parser extracts PLUGIN_NAME, VERSION, AUTHOR, DESCRIPTION._

- assert info['name'] == "My Cool Plugin"
- assert info['version'] == "1.2.3"
- assert info['author'] == "Test Author"
- assert info['description'] == "Does amazing things"
- assert info['category'] == "Analysis"

### test_plugin_info_version_tuple
_Version tuples like (1, 0, 0) should be joined as '1.0.0'._

- assert info['version'] == "2.3.1"

### test_plugin_info_fallback_dunder
___version__ and __author__ should be used when PLUGIN_ variants are missing._

- assert info['version'] == "0.5.0"
- assert info['author'] == "Dunder Author"

### test_plugin_info_docstring_fallback
_Module docstring should be used as description if PLUGIN_DESCRIPTION is missing._

- assert info['description'] == "This is the plugin description from docstring."

### test_plugin_info_missing_file
_Non-existent file should return defaults, not crash._

- assert info['name'] == "nonexistent.py"
- assert info['version'] == "Unknown"

### test_plugin_info_syntax_error
_File with syntax errors should return defaults, not crash._

- assert info['version'] == "Unknown"

### test_plugin_info_empty_file
_Empty file should return defaults._

- assert info['name'] == "empty.py"

### test_register_menu_action
_register_menu_action should store action metadata._

- assert len(pm.menu_actions) == 1
- assert pm.menu_actions[0]['plugin'] == "TestPlugin"

### test_register_toolbar_action
_No description_

- assert len(pm.toolbar_actions) == 1

### test_register_export_action
_No description_

- assert len(pm.export_actions) == 1
- assert pm.export_actions[0]['label'] == "Export as PDF"

### test_register_optimization_method
_No description_

- assert "MMFF94" in pm.optimization_methods

### test_register_file_opener
_No description_

- assert ".cif" in pm.file_openers

### test_register_file_opener_priority
_Higher priority opener should replace lower priority one in sorted order._

- assert pm.file_openers[".xyz"][0]['callback'] == cb_high

### test_register_analysis_tool
_No description_

- assert len(pm.analysis_tools) == 1

### test_register_save_handler
_No description_

- assert "TestPlugin" in pm.save_handlers

### test_register_load_handler
_No description_

- assert "TestPlugin" in pm.load_handlers

### test_register_3d_style
_No description_

- assert "Wireframe" in pm.custom_3d_styles

### test_register_document_reset_handler
_No description_

- assert len(pm.document_reset_handlers) == 1

### test_invoke_document_reset_handlers
_All registered reset handlers should be called._

- assert called == ["P1", "P2"]

### test_discover_plugins_empty_dir
_Empty plugin directory should return no plugins._

- assert plugins == []

### test_discover_plugins_single_file
_Single .py file in plugin dir should be discovered._

- assert len(plugins) >= 1
- assert "Hello" in names

### test_discover_plugins_package
_Folder with __init__.py should be discovered as package plugin._

- assert len(plugins) >= 1
- assert "Package Plugin" in names

### test_discover_plugins_ignores_dunder_files
_Files starting with __ should be ignored, and root __init__.py shouldn't count as a plugin._

- assert len(plugins) == 1 # The 'category' package itself
- assert "__trash" not in names

### test_ensure_plugin_dir_creates_directory
_ensure_plugin_dir should create the directory if it doesn't exist._

- assert os.path.isdir(new_dir)

## tests/integration/test_calculation_worker.py

### test_calculation_worker_bond_length_validation
_Integration test: Run a real 3D optimization for Ethane.
Compare the resulting C-C bond length against RDKit's reference calculated value.
No hard-coded bond length value is used._

- assert mol_3d is not None
- assert len(cc_bonds) >= 1
- assert dist_measured == pytest.approx(dist_ref, rel=1e-2)

### test_calculation_worker_angle_validation
_Integration test: Run 3D optimization for Propane.
Compare the C-C-C angle against RDKit's reference._

- assert central != -1
- assert angle_measured == pytest.approx(angle_ref, abs=5.0), \

### test_calculation_worker_dihedral_validation
_Integration test: Run 3D optimization for n-Butane.
Check if it converges to a staggered conformation (Trans or Gauche).
No hardcoded values; check against RDKit staggered preference._

- assert is_staggered

### test_calculation_worker_conversion_no_optimize
_Integration test: 2D to 3D without optimization.
Check if atom counts and symbols are preserved._

- assert "F" in symbols
- assert "Cl" in symbols
- assert mol_3d.GetNumAtoms() == 2 # F-Cl without Hs (default for this simple case if not specified)
- assert dist > 0.1

## tests/unit/test_app_logic.py

### test_worker_halt_logic
_Verify that CalculationWorker respects the halt flag (software flow)._

- assert "Halted" in str(blocker.args[0])

### test_ez_preservation_logic
_Verify worker preserves explicit E/Z labels even if RDKit might lose them._

- assert bond.GetStereo() == Chem.BondStereo.STEREOZ

### test_molecular_data_fallback_serialization
_Verify MolecularData uses manual string construction if RDKit fails (fallback logic)._

- assert mol_block is not None
- assert "MoleditPy" in mol_block
- assert "X  " in mol_block
- assert "C  " in mol_block

### test_coordinate_mapping_primary
_Verify primary RDKit conversion logic uses 'pos' attribute correctly._

- assert len(atom_lines) == 2
- assert abs(x_coord - expected_x) < 0.0001

### test_coordinate_mapping_fallback
_Verify fallback serialization (reverted logic) uses atom['item'].pos()._

- assert mol_block is not None
- assert "MoleditPy" in mol_block # Check header of fallback block
- assert len(atom_lines) == 2
- assert abs(x_coord - expected_x) < 0.0001

### test_worker_direct_mode
_Verify worker handles 'direct' mode (software branching)._

- assert mol.GetNumAtoms() == 5

## tests/unit/test_geometry.py

### test_3d_bond_lengths
_Verify optimized 3D coordinates yield physical bond lengths._

- assert mol is not None
- assert cc_bond is not None
- assert 1.50 < dist < 1.60

### test_mirror_dialog_logic
_Verify MirrorDialog correctly manipulates coordinates and UI (software logic)._

- assert x_after == pytest.approx(-x_before, abs=1e-3)
- assert main_window.draw_molecule_3d.called
- assert main_window.update_chiral_labels.called
- assert main_window.push_undo_state.called

### test_planarize_logic
_Verify planarize functionality using the actual PlanarizeDialog logic (SVD-based coplanarity check)._

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
_No description_

- assert os.path.exists(project_file)
- assert len(io_handler.data.atoms) == 1
- assert next(iter(io_handler.data.atoms.values()))['charge'] == 1

### test_xyz_export_logic
_No description_

- assert os.path.exists(xyz_file)
- assert xyz_mol.GetNumAtoms() == mol.GetNumAtoms()

## tests/unit/test_parsers.py

### test_fix_mol_block
_No description_

- assert "V2000" in lines[3]
- assert len(lines[3]) >= 39

### test_load_mol_file_logic
_No description_

- assert len(parser.data.atoms) == 2
- assert any(a['symbol'] == 'C' for a in parser.data.atoms.values())

### test_xyz_parsing_logic
_No description_

- assert mol is not None
- assert mol.GetNumAtoms() == 3
- assert any(a.GetSymbol() == 'O' for a in mol.GetAtoms())

### test_load_xyz_file_with_estimation
_No description_

- assert mol is not None
- assert mol.GetNumAtoms() == 2
- assert mol.GetNumBonds() == 1

### test_save_as_mol_logic
_No description_

- assert os.path.exists(save_path)
- assert "C  " in content

## tests/unit/test_properties.py

### test_analysis_window_regular_mol
_Verify AnalysisWindow calculates and displays properties for regular molecules._

- assert formula_found
- assert "C2H6O" in value

### test_analysis_window_xyz_derived
_Verify AnalysisWindow uses manual logic for XYZ-derived structures._

- assert "C2HO" in formula_val
- assert not smiles_present # SMILES should be withheld for XYZ

## tests/unit/test_stereochemistry.py

### test_wedge_dash_mapping
_Verify that Wedge/Dash in MolecularData mapped to RDKit BondDir._

- assert mol is not None
- assert bond.GetBondDir() == Chem.BondDir.BEGINWEDGE
- assert bond.GetBondDir() == Chem.BondDir.ENDWEDGE # If swapped
- assert bond.GetBondDir() == Chem.BondDir.BEGINDASH
- assert bond.GetBondDir() == Chem.BondDir.ENDDASH # If swapped

### test_ez_stereo_persistence
_Verify E/Z double bond configurations are maintained._

- assert bond.GetStereo() == Chem.BondStereo.STEREOE
- assert bond_z.GetStereo() == Chem.BondStereo.STEREOZ

### test_chiral_r_s_consistency
_Verify that R/S markers match RDKit descriptors for a known chiral center._

- assert atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED

### test_stereo_confirmation
_Verify that mirroring a molecule actually inverts its stereochemistry using MirrorDialog._

- assert mol.GetAtomWithIdx(1).GetProp("_CIPCode") == "R" # S becomes R

### test_stereo_loss_on_planarize
_Verify that planarizing a chiral center removes its chirality using PlanarizeDialog._

- assert not mol.GetAtomWithIdx(1).HasProp("_CIPCode")

