import os
import glob
import math
from pymol import cmd, cgo

def draw_cgo_arrow(obj_name, arrow_name):
    """
    Extracts the first and last coordinates of a loaded PyMOL object 
    and draws a 3D CGO arrow between them.
    """
    try:
        # Get coordinates of all atoms in the trajectory
        coords = cmd.get_coords(obj_name)
        if coords is None or len(coords) < 2:
            return
            
        p1 = coords[0]  # Entry point
        p2 = coords[-1] # Exit point
        
        # Calculate distance and direction
        dx, dy, dz = p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]
        length = math.sqrt(dx**2 + dy**2 + dz**2)
        if length == 0: 
            return
            
        direction = (dx/length, dy/length, dz/length)
        cone_length = min(2.0, length * 0.2) 
        
        cyl_end = (
            p2[0] - direction[0] * cone_length,
            p2[1] - direction[1] * cone_length,
            p2[2] - direction[2] * cone_length
        )
        
        # Arrow aesthetics
        radius = 0.3
        cone_radius = 0.8
        r, g, b = 1.0, 0.8, 0.2  # Gold/Orange color for high visibility
        
        # Build the CGO array
        arrow_cgo = [
            cgo.CYLINDER, 
            p1[0], p1[1], p1[2], 
            cyl_end[0], cyl_end[1], cyl_end[2], 
            radius, 
            r, g, b, r, g, b,
            
            cgo.CONE, 
            cyl_end[0], cyl_end[1], cyl_end[2], 
            p2[0], p2[1], p2[2], 
            cone_radius, 0.0, 
            r, g, b, r, g, b, 
            1, 1 
        ]
        
        cmd.load_cgo(arrow_cgo, arrow_name)
    except Exception as e:
        print(f"Warning: Could not draw arrow for {obj_name}. Error: {e}")

def visualize_kinetics():
    """
    Main function to load the crystal structure, apply formatting, 
    and loop through all water trajectory PDBs.
    """
    print("Initializing Kinetic Visualization...")

    # --- 1. Load and Format the Base Structure ---
    base_pdb = "/Users/kanishkkondaka/Desktop/margaret_files/start_Cl_sim_new_unit_cell.pdb"
    base_name = "start_Cl_sim_new_unit_cell"
    
    # Load reference structure
    cmd.load(base_pdb, base_name)
    
    # Setup selections and display settings
    cmd.select(f"{base_name}_no_SOL", f"{base_name} and not resn SOL")
    cmd.hide("everything", f"{base_name}_no_SOL")
    
    cmd.show("spheres", "resn OEX")
    cmd.label("resn OEX", "name")
    cmd.set("label_size", 9)

    # --- 2. Process the Water Trajectories ---
    water_dir = "/Users/kanishkkondaka/Desktop/margaret_files/waters_of_interest"
    
    # Find all PDB files in the target directory
    pdb_files = glob.glob(os.path.join(water_dir, "*.pdb"))
    
    if not pdb_files:
        print(f"Error: No PDB files found in {water_dir}")
        return
        
    print(f"Found {len(pdb_files)} water trajectories. Loading...")

    for pdb_path in pdb_files:
        # Extract the filename without the .pdb extension to use as the object name
        obj_name = os.path.basename(pdb_path).replace(".pdb", "")
        
        # Load the water trajectory
        cmd.load(pdb_path, obj_name)
        
        # Apply coloring based on atomic index (rank)
        cmd.spectrum("rank", "blue_white_red", obj_name)
        
        # Draw the custom CGO arrow
        arrow_name = f"arrow_{obj_name}"
        draw_cgo_arrow(obj_name, arrow_name)
        
        # Optional: display the waters as spheres for better visibility
        cmd.show("spheres", obj_name)
        cmd.set("sphere_scale", 0.3, obj_name)

    print("Visualization complete!")

# Automatically run the main function when the script is executed in PyMOL
visualize_kinetics()
