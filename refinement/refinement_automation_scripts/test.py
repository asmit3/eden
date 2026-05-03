import sys
import os
import numpy as np
from Bio.PDB import PDBParser
from refine1_script import run_main

def get_latest_pdb_file(directory):
    pdb_files = [f for f in os.listdir(directory) if f.endswith(".pdb")]
    pdb_files.sort(key=lambda f: os.path.getmtime(os.path.join(directory, f)))
    return os.path.join(directory, pdb_files[-1])  

def get_first_atom_xyz(pdb_file): #change to cctbx
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("refined", pdb_file)
    first_atom = next(structure.get_atoms())
    return first_atom.coord

def main():
    run_main()

    output_dir = "/Users/lbnl/Desktop/Srikala_Munukutla_DesktopFiles/refinement_1/script_run"
    latest_pdb = get_latest_pdb_file(output_dir)
    expected_filename = "refine_001.pdb"

    assert os.path.basename(latest_pdb) == expected_filename

    expected_xyz = np.array([-9.079, 4.689, 5.833])
    actual_xyz = get_first_atom_xyz(latest_pdb)

    assert np.allclose(actual_xyz, expected_xyz, atol=1e-3)

    """for filename, expected_xyz in zip(expected_files, expected_coords):
        pdb_path = get_pdb_path(output_dir, filename)
        actual_xyz = get_first_atom_xyz(pdb_path)
        assert np.allclose(actual_xyz, expected_xyz, atol=1e-3)"""

if __name__ == "__main__":
    main()







