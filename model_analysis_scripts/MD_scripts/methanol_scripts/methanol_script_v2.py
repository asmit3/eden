import numpy as np
import random
import argparse
import iotbx.pdb
from scitbx.array_family import flex

VDW_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, 
    "S": 1.80, "P": 1.80, "F": 1.47, "CL": 1.75
}

def get_vdw(element):
    element = element.strip().upper()
    return VDW_RADII.get(element, 1.5)

def check_clash(query_sites, query_radii, obstacle_coords, obstacle_radii):
    if len(obstacle_coords) == 0:
        return False
        
    for i, q_site in enumerate(query_sites):
        q_rad = query_radii[i]
        diff = obstacle_coords - q_site
        dists_sq = diff.dot() 
        
        thresholds = 0.8 * (obstacle_radii + q_rad) # using square distances for faster computation, square roots take a lot of time to compute
        thresholds_sq = thresholds * thresholds
        
        if (dists_sq < thresholds_sq).count(True) > 0:
            return True
    return False


def extract_water_xyz(input_pdb):
    pdb_inp = iotbx.pdb.input(file_name = input_pdb)
    hierarchy = pdb_inp.construct_hierarchy() 

    selection_str = "(resname SOL or resname WAT or resname HOH) and (name O or name OW)"

    cache = hierarchy.atom_selection_cache()
    selection = cache.selection(selection_str)

    water_oxygen_coords = hierarchy.atoms().select(selection).extract_xyz()

    coords_list = list(water_oxygen_coords)
    return coords_list


def generate_methanols(seed, pdb_inp, output_file):
    random.seed(seed)

    symmetry = pdb_inp.crystal_symmetry()
    uc_params = symmetry.unit_cell().parameters()
    box_x = uc_params[0]
    box_y = uc_params[1]
    box_z = uc_params[2]

    num_methanols = 2322

    sys_hierarchy = pdb_inp.construct_hierarchy()# hierarchy
    
    cache = sys_hierarchy.atom_selection_cache()
    prot_sel = cache.selection("not (resname SOL or resname WAT or resname HOH)") # selection syntax
    prot_atoms = sys_hierarchy.atoms().select(prot_sel)
    
    obs_coords = prot_atoms.extract_xyz()
    obs_radii = flex.double([get_vdw(a.element) for a in prot_atoms])

    methanol_pdb = iotbx.pdb.input(file_name="/Users/yyklab/Downloads/methanol.pdb")
    base_hierarchy = methanol_pdb.construct_hierarchy()
    base_atoms = base_hierarchy.atoms()
    base_sites = base_atoms.extract_xyz()
    meth_radii = [get_vdw(a.element) for a in base_atoms]

    out_hierarchy = sys_hierarchy.deep_copy()
    out_model = out_hierarchy.models()[0] 

    # we will store just the methanol coordinates here to check against waters later
    meth_only_coords = flex.vec3_double()
    meth_only_radii = flex.double()

    placed = 0
    attempts = 0
    max_attempts = 1000000 

    print(f"attempting to place {num_methanols} methanols...")

    while placed < num_methanols and attempts < max_attempts:
        attempts += 1
        
        dx = random.uniform(0, box_x)
        dy = random.uniform(0, box_y)
        dz = random.uniform(0, box_z)
        
        translated_sites = base_sites + (dx, dy, dz)
        
        # check against protein AND previously placed methanols
        if not check_clash(translated_sites, meth_radii, obs_coords, obs_radii):
            methanol_copy = base_hierarchy.deep_copy()
            methanol_copy.atoms().set_xyz(translated_sites)
            
            for chain in methanol_copy.models()[0].chains():
                chain_copy = chain.detached_copy()
                for rg in chain_copy.residue_groups():
                    rg.resseq = str(placed + 1)
                out_model.append_chain(chain_copy)

            # add to the general obstacle pool (so future methanols don't clash)
            obs_coords.extend(translated_sites)
            obs_radii.extend(flex.double(meth_radii))
            
            # add to the methanol-only pool (so we can check waters later)
            meth_only_coords.extend(translated_sites)
            meth_only_radii.extend(flex.double(meth_radii))
            
            placed += 1
            if placed % 500 == 0:
                print(f"placed {placed} / {num_methanols} methanols...")

    if placed < num_methanols:
        print(f"warning: only placed {placed} methanols due to space constraints.")
    else:
        print("success! all methanols placed.")

    print("scanning for overlapping water molecules (check oxygens only)...")
    waters_removed = 0
    
    for model in out_hierarchy.models():
        for chain in model.chains():
            rgs_to_remove = []
            
            for rg in chain.residue_groups():
                resname = rg.atom_groups()[0].resname.strip()
                if resname in ["SOL", "WAT", "HOH"]:
                    
                    # create flex arrays just for the oxygen(s) in this water molecule
                    ox_sites = flex.vec3_double()
                    ox_radii = flex.double()
                    
                    for atom in rg.atoms():
                        # identify oxygen by element symbol or common atom names
                        if atom.element.strip().upper() == "O" or atom.name.strip() in ["O", "OW", "OH2"]:
                            ox_sites.append(atom.xyz)
                            ox_radii.append(get_vdw(atom.element))
                    
                    # if we found an oxygen and it clashes
                    if len(ox_sites) > 0 and check_clash(ox_sites, ox_radii, meth_only_coords, meth_only_radii):
                        rgs_to_remove.append(rg)
            
            for rg in rgs_to_remove:
                chain.remove_residue_group(rg)
                waters_removed += 1

    print(f"cleanup complete: removed {waters_removed} clashing water molecules.")

    out_hierarchy.write_pdb_file(file_name=output_file, crystal_symmetry=symmetry)


if __name__ == "__main__": 
    parser = argparse.ArgumentParser(description="Add random methanols to a system and remove clashing waters.")
    parser.add_argument("-f", "--file", required=True, help="Input PDB file (Receptor/Protein + Waters)")
    parser.add_argument("-o", "--out", required=True, help="Output PDB file")
    parser.add_argument("-s", "--seed", type=int, default=42, help="Random seed for reproducibility")
    
    args = parser.parse_args()
    
    pdb_input_obj = iotbx.pdb.input(file_name=args.file)
    generate_methanols(args.seed, pdb_input_obj, args.out)