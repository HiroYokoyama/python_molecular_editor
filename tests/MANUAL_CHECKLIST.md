# Software Verification Checklist

This comprehensive checklist is designed to verify that **MoleditPy** functions as documented, matches the intended design, and maintains scientific integrity.

## 1. Installation & Launch
- [ ] **Package Installation**: Does the installation command (e.g., `pip install .`) execute successfully?
- [ ] **Shortcut Creation**: Does the post-installation script (`moleditpy-installer`) correctly create a shortcut in the application menu (Start Menu/Applications)?
- [ ] **CLI Launch**: Does the primary launch command open the application from the terminal?
- [ ] **Standalone Installer**: Does the standalone installer (exe) work without a pre-existing environment?
- [ ] **Dependency Check**: Does the application launch without errors regarding core scientific or geometry libraries?

## 2. Interface Layout
- [ ] **Main Components**: Are the Menu Bar, Main Toolbar, Template Toolbar, 2D Edit View, 3D View, and Status Bar all visible?
- [ ] **Splitter Functionality**: Can the boundary between the 2D and 3D views be resized by dragging?
- [ ] **Layout Shortcuts**: Do the layout presets (e.g., `Ctrl+1`, `Ctrl+2`, `Ctrl+3`) resize the panels to the correct ratios?

## 3. 2D Editing Functions

### Basic Drawing
- [ ] **Atom Placement**: Can atoms be placed using toolbar buttons or keyboard shortcuts?
- [ ] **Bond Creation**: Can bonds be created by dragging from one atom to another?
- [ ] **Sprouting**: Can a new atom be created by dragging from an existing atom into empty space?
- [ ] **Bond Modification**: Does clicking an existing bond in single/double/triple bond mode change bond order?
- [ ] **Element Switching**: Can an existing atom be changed by selecting a new element or using hover-key shortcuts?

### Properties & Stereochemistry
- [ ] **Charge Adjustment**: Do the +/- buttons adjust the formal charge, and does a right-click reset it?
- [ ] **Radical State**: Does the Radical button toggle the number of electrons (0 -> 1 -> 2 -> 0)?
- [ ] **Stereo Bonds**: Can single bonds be changed to Wedge or Dash representations?
- [ ] **E/Z Configuration**: Can double bonds be toggled between E/Z stereo configurations?

### Editing Tools
- [ ] **Selection**: Can items be selected individually, via Shift-click, and via rectangular drag?
- [ ] **Deletion**: Does the right-click action (or Delete/Backspace key) remove selected items?
- [ ] **Undo/Redo**: Do the Undo/Redo functions work correctly for all 2D drawing operations?
- [ ] **Copy/Paste**: Do the Copy/Paste functions work within the canvas?
- [ ] **Clean Up**: Does the 2D Clean Up function correctly organize the structure coordinates?
- [ ] **Hydrogen Management**: Do the "Add Hydrogens" and "Remove Hydrogens" functions update the graph correctly?
- [ ] **Clear All**: Does the "Clear All" function successfully wipe the canvas?

### Templates
- [ ] **Standard Templates**: Can standard rings (e.g., Benzene) be placed and merged with existing structures?
- [ ] **User Templates**: Can the current structure be saved as a custom user template?
- [ ] **Retrieval**: Can saved user templates be retrieved and placed on the canvas?

## 4. File Operations

### Project Management
- [ ] **Save Project**: Can the workspace be saved as a project file (`.pmeprj`) with all data intact?
- [ ] **Open Project**: Can a saved project file be re-opened with both 2D and 3D states preserved?

### Import
- [ ] **Chemical Formats**: Does importing standard files (e.g., MOL, SDF) load the structure into the 2D view correctly?
- [ ] **SMILES**: Does entering a SMILES string generate the expected 2D topology?
- [ ] **InChI**: Does entering an InChI string generate the expected 2D topology?
- [ ] **3D Direct Import**: Do specific 3D formats load directly into the 3D view, bypassing the 2D editor?

### Export
- [ ] **Images**: Can the structure be exported as PNG (2D/3D) and SVG (2D) images?
- [ ] **Chemical Data**: Can the data be exported into valid MOL and XYZ files?
- [ ] **3D Models**: Can the 3D model be exported as an STL file (suitable for 3D printing)?
- [ ] **Color Preservation**: Does the OBJ/MTL export preserve atom color information?

## 5. 3D Generation & Visualization
- [ ] **2D to 3D Conversion**: Does the conversion function generate a valid 3D structure from a 2D drawing?
- [ ] **Structural Optimization**: Does the 3D optimization function refine the geometry using the selected force field?
- [ ] **View Controls**: Do Zoom (Wheel), Rotate (Drag), and Pan (Shift+Drag) operate smoothly?
- [ ] **Display Styles**: Can the user switch between Ball & Stick, CPK, Wireframe, and Stick styles?
- [ ] **Aromatic Representation**: Can aromatic rings be toggled between Kekulé (double bonds) and Torus (circles)?

## 6. 3D Measurement & Editing

### Measurement
- [ ] **Distance**: Selecting 2 atoms displays the correct interatomic distance.
- [ ] **Angle**: Selecting 3 atoms displays the correct bond angle.
- [ ] **Dihedral**: Selecting 4 atoms displays the correct dihedral angle.

### 3D Editing (Manual)
- [ ] **Atom Drag**: Can atoms be moved manually in 3D space using the designated modifier key?
- [ ] **Global Translation**: Does the Translation tool move the entire molecule correctly?
- [ ] **Axis Alignment**: Can the molecule be aligned to X/Y/Z axes based on selected atoms?
- [ ] **Plane Alignment**: Can the molecule be aligned to a specific plane based on a selection of 3+ atoms?
- [ ] **Mirror Image**: Can a mirror image of the molecule be generated across a specified plane?
- [ ] **Precision Adjustments**: Can bond lengths, angles, and dihedrals be set to specific numeric values?
- [ ] **Planarize**: Does the function flatten selected atoms onto a best-fit plane?

### Constrained Optimization
- [ ] **Constraint Creation**: Can Distance/Angle/Dihedral constraints be added and listed?
- [ ] **Visual Feedback**: Does selecting a constraint highlight the corresponding atoms in the 3D view?
- [ ] **Execution**: Does the optimization respect the defined constraints during the calculation?

## 7. Analysis & Information
- [ ] **Atom Labels**: Can labels (Index, Coordinates, Element) be toggled on/off in the 3D view?
- [ ] **Chiral Labels**: Does the "Show Chiral Labels" option correctly display R/S markers in 3D viewer?
- [ ] **Analysis Properties**: Does the analysis dialog display correct properties (Molecular Weight, Formula, LogP, TPSA)?
- [ ] **Data Access**: Can individual property values be copied to the clipboard?

## 8. Settings & Customization
- [ ] **Appearance**: Can the background and bond colors be customized successfully?
- [ ] **3D Rendering**: Can lighting intensity and material properties (e.g., shininess) be adjusted?
- [ ] **Global Reset**: Does the "Reset All Settings" function restore the application to its default state?

## 9. Verification & Metadata
- [ ] **Version Info**: Does the "About" dialog display the correct current version?
- [ ] **Licensing**: Is the license information correctly listed and accessible?
- [ ] **Project Integrity**: Does the provided repository link point to the correct source?

---

### Verification Sign-off

| Field | Value |
| :--- | :--- |
| **Version Verified** | - |
| **Date** | - |
| **Lead Developer** | - |
