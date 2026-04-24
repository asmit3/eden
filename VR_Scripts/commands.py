import os
from chimerax.core.commands import run as chimerax_run

def set(session):
    script = os.path.dirname(os.path.realpath(__file__))
    chimerax_run(session, "close all")
    chimerax_run(session, "set bgColor white")
    chimerax_run(session, f"open '{os.path.join(script, '0F.pdb')}'") 
    chimerax_run(session, f"open '{os.path.join(script, '0F_2mFo-DFc.ccp4')}'") 
    chimerax_run(session, "hide all")
    chimerax_run(session, "volume #2 hide")
    chimerax_run(session, "show /A:170,189,601 atoms") #targets
    chimerax_run(session, "style /A:170,189,601 stick")

    chimerax_run(session, "volume zone #2 nearAtoms /A:170,189 range 1.5 newMap true") #tight zone, reduced from 2.5 to 1.5
    
    #turns density into blobs/spheres
    chimerax_run(session, "volume #3 style surface")
    chimerax_run(session, "volume #3 level 4.5") #higher number = tighter blobs
    chimerax_run(session, "color #3 skyblue")
    chimerax_run(session, "volume #3 transparency 0.2") 
    
    chimerax_run(session, "view /A:601")
    chimerax_run(session, "zoom 2.0")

set(session)