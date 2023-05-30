Scripts:

H_bond_finder_final

- Given any atom in a PDB, automatically plot the potential surrounding hydrogen bonding environment
- Only takes one PDB file at once, but generates thorough depiction enbironment with all bonding atoms and parent residues
- Returns dataframe with distances in text format for given file and center atom listing interactions

H_bond_finder_minimalistic

- Given a list of atoms, automatically plot a truncated version of the hydrogen bonding environment
- Capabale of analysing the bonding networks of multiple atoms across multiple timepoints (multiple PDBs)
- Returns dataframe with distances for bonding of target atom across all provided PDBs
- Minimalistic plots used to better study geometry – ideal tetrahedron for comparison toggled on and off with "t" key

Future directions for the H-bond network finder project:

- incorporate analysis of bond angles 
- compare geometry of bonding interactions to ideal tetrahedron 
- make movies with scenes files over different timepoints
- programmable selection of viewing angle
- automatically identify patterns in H-bond lengths, angles, and presences?

Potential features to fix:

- verify if conformation filtering is working
- figure out what to do about OEC vs OEZ inclusion (maybe this can just be left to the user?)
