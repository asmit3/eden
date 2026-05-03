import os
from phenix.command_line import superpose_pdbs
from phenix.programs import phenix_refine
from phenix.programs.phenix_refine import custom_process_arguments
from iotbx.data_manager import DataManager
from iotbx.cli_parser import run_program

# --- CONFIGURATION ---
step1_pdb = "step1.pdb"  # Already refined PDB from Step 1
template_pdb = "template.pdb"  # Template file to align
atoms_to_add = "metals"  # Choose: "all", "metals", "oxygen"

# --- Step 2a: Align Template to Step1 ---
superposed_pdb = "superposed.pdb"
var1 = f"input.pdb_file_name_fixed={step1_pdb}"
var2 = f"input.pdb_file_name_moving={template_pdb}"
var3 = f"output.file_name={superposed_pdb}"
superpose_pdbs.run([var1, var2, var3])

# --- Step 2b: Add OEC Atoms ---
step2_input_pdb = "step2_input.pdb"
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
    selected_elements = allowed_elements[element_selection]
    
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


add_OEC_atoms(step1_pdb, superposed_pdb, step2_input_pdb, element_selection=atoms_to_add)

# --- Step 2c: Refine After OEC Addition ---
def run_refine(input_pdb, prefix):
    args = [
        input_pdb,
        f"output.prefix={prefix}",
        "refinement.main.number_of_macro_cycles=3"
    ]
    run_program(
        program_class=phenix_refine.Program,
        args=args,
        custom_process_arguments=custom_process_arguments,
        logger=None
    )

run_refine(step2_input_pdb, prefix="step2_refine")

# --- Step 3: Add Solvent and Refine ---
step3_input_pdb = "step3_input.pdb"
def add_solvent(input_pdb, output_pdb):
    args = [
        input_pdb,
        f"output.prefix={os.path.splitext(output_pdb)[0]}",
        "refinement.main.number_of_macro_cycles=1",
        "refinement.water=true"
    ]
    run_program(
        program_class=phenix_refine.Program,
        args=args,
        custom_process_arguments=custom_process_arguments,
        logger=None
    )

add_solvent("step2_refine_001.pdb", step3_input_pdb)

# --- Final Refinement After Solvent Addition ---
run_refine(step3_input_pdb, prefix="step3_refine")