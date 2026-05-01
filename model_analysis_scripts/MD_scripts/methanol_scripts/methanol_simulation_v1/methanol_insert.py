import random
import argparse
import iotbx.pdb

def remove_random_waters(input_pdb, output_pdb, num_to_remove, seed):
    pdb_inp = iotbx.pdb.input(file_name=input_pdb)
    hierarchy = pdb_inp.construct_hierarchy()


    water_residue_groups = []

# construct_hierarchy() creates this: Model → Chain → Residue Group → Atom Group → Atom

    for model in hierarchy.models():
        for chain in model.chains():
            for rg in chain.residue_groups(): # residue groups
                resname = rg.atom_groups()[0].resname.strip()
                if resname in ["SOL", "WAT", "HOH"]:
                    water_residue_groups.append((chain, rg))

    
    total_waters = len(water_residue_groups)
    print(f"there are {total_waters} waters")

    if total_waters < num_to_remove:
        print(f"too many waters, {num_to_remove} is more than the total number ofw waters: {total_waters}.")
        return

    random.seed(seed)
    waters_to_remove = random.sample(water_residue_groups, num_to_remove)

    for chain, rg in waters_to_remove:
        chain.remove_residue_group(rg)


    with open(output_pdb, "w") as f_out:
        f_out.write(hierarchy.as_pdb_string(crystal_symmetry=pdb_inp.crystal_symmetry()))

    print(f"removed {num_to_remove} waters")
    print(f"saved to {output_pdb}")


if __name__ == "__main__": # https://docs.python.org/3/howto/argparse.html
    parser = argparse.ArgumentParser(description="Remove random waters using iotbx.")
    parser.add_argument("-f", "--file", required=True, help="Input PDB file")
    parser.add_argument("-o", "--out", required=True, help="Output PDB file")
    parser.add_argument("-n", "--num", type=int, default=60, help="Number of waters to delete")
    parser.add_argument("-s", "--seed", type=int, default=42, help="Random seed for reproducibility")
    
    args = parser.parse_args()
    remove_random_waters(args.file, args.out, args.num, args.seed)