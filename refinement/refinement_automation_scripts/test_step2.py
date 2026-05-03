import os
from phenix.command_line import superpose_pdbs
from iotbx.data_manager import DataManager
from IPython import embed

def add_OEC_atoms(step1_file, aligned_template_file, output_pdb, element_selection="all"):
    dm = DataManager()
    dm.set_overwrite(True)
    model_step1 = dm.get_model(step1_file)
    model_template = dm.get_model(aligned_template_file)
    ph1 = model_step1.get_hierarchy()
    ph2 = model_template.get_hierarchy()

    allowed_elements = {
        "all": {"MN", "CA", "O"},
        "metals": {"MN", "CA"},
        "oxygen": {"O"}
    }
    selected_elements = allowed_elements.get(element_selection.lower())

    for chain in ph2.chains():
        if chain.id.strip() in ["A", "a"]:
            for rg in chain.residue_groups():
                if rg.resid().strip() == "601":
                    rg_copy = rg.detached_copy()
                    for atom_group in rg_copy.atom_groups():
                        delete_atoms = [atom for atom in atom_group.atoms() if atom.element.upper() not in selected_elements]
                        for atom in delete_atoms:
                            atom_group.remove_atom(atom)
                    for target_chain in ph1.chains():
                        if target_chain.id == chain.id:
                            target_chain.append_residue_group(rg_copy)
                            break

    dm.write_model_file(model_step1, output_pdb)

step1_output = os.path.abspath("step1.pdb")
template_pdb = os.path.abspath("template.pdb")
superposed_path = os.path.abspath("superposed.pdb")
step2_input = os.path.abspath("step2_input.pdb")

var1 = f"input.pdb_file_name_fixed = {step1_output}"
var2 = f"input.pdb_file_name_moving = {template_pdb}"
var3 = f"output.file_name = {superposed_path}"
print(var1)
print(var2)
print(var3)
result = superpose_pdbs.run([var1, var2, var3])

choice = input("Enter choice [1/2/3]: ")

if choice == "1":
    atoms_to_add = "all"
elif choice == "2":
    atoms_to_add = "metals"
elif choice == "3":
    atoms_to_add = "oxygen"

add_OEC_atoms(step1_output, superposed_path, step2_input, element_selection=atoms_to_add)
