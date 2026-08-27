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
        coords = cmd.get_coords(obj_name)
        if coords is None or len(coords) < 2:
            return
            
        p1 = coords[0]  
        p2 = coords[-1] 
        
        dx, dy, dz = p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]
        length = math.sqrt(dx**2 + dy**2 + dz**2)
        if length == 0: 
            return
            
        direction = (dx/length, dy/length, dz/length)
        cone_length = min(1.5, length * 0.2) 
        
        cyl_end = (
            p2[0] - direction[0] * cone_length,
            p2[1] - direction[1] * cone_length,
            p2[2] - direction[2] * cone_length
        )
        
        # --- SMALLER ARROW AESTHETICS ---
        radius = 0.5        # Shaft thickness (reduced from 0.3)
        cone_radius = 1.0    # Head thickness (reduced from 0.8)
        r, g, b = 1.0, 0.8, 0.2  
        
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
    print("Initializing Kinetic Visualization...")

    # --- 1. Load and Format the Base Structure ---
    base_pdb = "/Users/kanishkkondaka/Desktop/margaret_files/start_Cl_sim_new_unit_cell.pdb"
    base_name = "start_Cl_sim_new_unit_cell"
    
    cmd.load(base_pdb, base_name)
    
    cmd.select(f"{base_name}_no_SOL", f"{base_name} and not resn SOL")
    cmd.hide("everything", f"{base_name}_no_SOL")
    
    cmd.show("spheres", "resn OEX")
    cmd.label("resn OEX", "name")
    cmd.set("label_size", 9)

# Set background to white so the colored waters pop
    cmd.bg_color("white")

    # --- 2. Process the Water Trajectories ---
    water_dir = "/Users/kanishkkondaka/Desktop/margaret_files/waters_of_interest"
    pdb_files = glob.glob(os.path.join(water_dir, "*.pdb"))
    
    if not pdb_files:
        print(f"Error: No PDB files found in {water_dir}")
        return
        
    print(f"Found {len(pdb_files)} water trajectories. Loading...")

    loaded_objects = []

    for pdb_path in pdb_files:
        obj_name = os.path.basename(pdb_path).replace(".pdb", "")
        arrow_name = f"arrow_{obj_name}"
        
        cmd.load(pdb_path, obj_name)
        
        # Force the representation to be ONLY spheres (clears out default lines/crosses)
        cmd.show_as("spheres", obj_name)
        
        # Set to full size (1.0) or slightly smaller (0.8) so they are highly visible
        cmd.set("sphere_scale", 0.5, obj_name)
        
        # Color chronologically
        cmd.spectrum("rank", "blue_white_red", obj_name)
        
        # Draw the arrow
        draw_cgo_arrow(obj_name, arrow_name)
        
        loaded_objects.append({'water': obj_name, 'arrow': arrow_name})
        
        # Disable them immediately so the screen isn't cluttered
        cmd.disable(obj_name)
        cmd.disable(arrow_name)

# --- 3. Automated Screenshot Loop ---
    image_dir = "/Users/kanishkkondaka/Desktop/margaret_files/water_images"
    os.makedirs(image_dir, exist_ok=True)
    
    print(f"Taking isolated screenshots and saving to: {image_dir}")
    
    # Ensure EVERYTHING else (like the protein base structure) is turned off
    cmd.disable("all")
    
    for obj in loaded_objects:
        w_name = obj['water']
        a_name = obj['arrow']
        
        # 1. Enable and explicitly show ONLY this water and arrow
        cmd.enable(w_name)
        cmd.enable(a_name)
        cmd.show("spheres", w_name)
        # 2. Zoom perfectly to the trajectory
        cmd.zoom(w_name)
        
        # 3. CRITICAL: Force PyMOL to update the graphics buffer before the photo
        cmd.refresh()
        cmd.enable(a_name)

        # 4. Take the screenshot
        img_path = os.path.join(image_dir, f"{w_name}.png")
        cmd.png(img_path, width=1920, height=1080, ray=0)
        
        # 5. Deselect and disable before the next water
        cmd.disable(w_name)
        cmd.disable(a_name)
        cmd.hide("spheres", w_name)

    print("All screenshots generated successfully!")

# Execute
visualize_kinetics()
