import os
from chimerax.core.commands import run as chimerax_run

def set(session):
    script = os.path.dirname(os.path.realpath(__file__))
    chimerax_run(session, "close all")
    chimerax_run(session, "set bgColor white")
    
    datasets = ["0F", "1F"]
    target_residue = "/A:601"

    for i, prefix in enumerate(datasets): #we have counter and the prefix
        base_id = (i * 3) + 1 #pdb file is #1, fo-fc map is #2, 2fo-fc map is #3 for OF files
        pdb_id = base_id
        m2fofc_id = base_id + 1
        mfofc_id = base_id + 2

        chimerax_run(session, f"open '{os.path.join(script, prefix + '.pdb')}'")
        chimerax_run(session, f"open '{os.path.join(script, prefix + '_2mFo-DFc.ccp4')}'")
        chimerax_run(session, f"open '{os.path.join(script, prefix + '_mFo-DFc.ccp4')}'")

        chimerax_run(session, f"hide all")
        chimerax_run(session, f"show #{pdb_id}/A,B,C,D cartoon")
        chimerax_run(session, f"show #{pdb_id}{target_residue} :<10 atoms") #select site surroundings and show as sticks
        chimerax_run(session, f"style #{pdb_id}{target_residue} :<10 stick")
        chimerax_run(session, f"show #{pdb_id}{target_residue} :<10 & :HOH atoms")
        chimerax_run(session, f"style #{pdb_id}{target_residue} :<10 & :HOH sphere") #show waters as spheres
        chimerax_run(session, f"hide #{pdb_id} atoms")
        chimerax_run(session, f"hide #{pdb_id} ribbons")
        chimerax_run(session, f"bond #{pdb_id}{target_residue}@MN4 #{pdb_id}{target_residue}@O5")
        chimerax_run(session, f"show #{pdb_id}{target_residue} atoms")
        chimerax_run(session, f"style #{pdb_id}{target_residue}@MN4,O5 ball")
        chimerax_run(session, f"size #{pdb_id}{target_residue}@MN4 atomRadius 1.2")
        chimerax_run(session, f"size #{pdb_id}{target_residue} stickRadius 0.1")
        chimerax_run(session, f"hide #{m2fofc_id},{mfofc_id} models")
        zone_id = f"{pdb_id}.1"
       

    if len(datasets) > 1:
        chimerax_run(session, "matchmaker #4 to #1")
    
    chimerax_run(session, f"view {target_residue}")
    chimerax_run(session, "zoom 1.2")

set(session)
