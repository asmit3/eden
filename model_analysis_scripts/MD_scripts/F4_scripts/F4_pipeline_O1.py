from pymol import cmd, cgo, CmdException
from pathlib import Path
import numpy as np
#----INPUTS-----
protein_file = "/Users/yyklab/Desktop/Lab_Files/restart_09_09_25/base_files/0F.pdb_fittedto55all.pdb" # 0F.pdb
H1_peaks = "/Users/yyklab/Desktop/Lab_Files/low_restraints/maps_peaks/peaks_low_restraints_justH1_1_rms_2.13.pdb"
H2_peaks = "/Users/yyklab/Desktop/Lab_Files/low_restraints/maps_peaks/peaks_low_restraints_justH2_1_rms_2.11.pdb"
OW_peaks = "/Users/yyklab/Desktop/Lab_Files/low_restraints/maps_peaks/peaks_low_restraints_justOW_1_rms_1.58.pdb"
water_ccp4_map = "/Users/yyklab/Desktop/Lab_Files/low_restraints/base_files/FFT_6/fcalc_lower_restraints_last55ns_p212121_justwater_1.ccp4"

#-------WATERS ARRAY--------
ticklabels2 = ['48.0', '21.0', '62.0', "'OEC_A_W1_601'", '62.0', "'ASN_A_OD1_181'", '52.0', '34.0', '26.0', '19.0', '21.0', '22.0', '28.0', "'ARG_A_NH1_334'", '66.0', "'ARG_A_NH1_334'", '66.0', "'OEC_A_W4_601'", '28.0', '30.0', '29.0', '29.0', '31.0', '33.0', '34.0', '37.0', '38.0', '32.0', '77.0', '38.0', '104.0']
ticklabels1 = ['ASP_A_O_61','ASP_A_OD2_61','SER_A_O_169','SER_A_O_169','GLY_A_O_171','ASN_A_ND2_181','ASN_A_O_335','PRO_A_O_340','OEC_A_O1_601','OEC_A_O4_601','OEC_A_W1_601','OEC_A_W2_601','OEC_A_W4_601','GLU_D_OE1_312','GLU_D_OE1_312','GLU_D_OE2_312','GLU_D_OE2_312',26.0,26.0,27.0,28.0,30.0,32.0,32.0,33.0,36.0,37.0,39.0,76.0,77.0,103.0]

ticklabels = ticklabels2 + ticklabels1 

def separate_array(input_array):
    integers_array = []
    strings_array = []
    
    for item in input_array:
        try:
            # First convert to float to handle the '.0', then convert to an integer
            num = int(float(item))
            integers_array.append(num)
        except ValueError:
            # If the conversion fails, it's a standard string (like your residue names)
            strings_array.append(item)
            
    return integers_array, strings_array


waters, res = separate_array(ticklabels)
res = [i.strip("'") for i in res]
print(res)
        

cmd.set_color("red1", [1.0,0.0,0.0])
cmd.set_color("red2", [1.0, 0.3, 0.3])
cmd.set_color("red3", [1.0,1.0,1.0])

def cgo_arrow(atom1='pk1', atom2='pk2', radius=0.5, gap=0.0, hlength=-1, hradius=-1, color='blue red', name=''):
    '''
DESCRIPTION

    Create a CGO arrow between two picked atoms.

ARGUMENTS

    atom1 = string: single atom selection or list of 3 floats {default: pk1}

    atom2 = string: single atom selection or list of 3 floats {default: pk2}

    radius = float: arrow radius {default: 0.5}

    gap = float: gap between arrow tips and the two atoms {default: 0.0}

    hlength = float: length of head

    hradius = float: radius of head

    color = string: one or two color names {default: blue red}

    name = string: name of CGO object
    '''

    from chempy import cpv

    
    radius, gap = float(radius), float(gap)
    hlength, hradius = float(hlength), float(hradius)

    try:
        color1, color2 = color.split()
    except:
        color1 = color2 = color
    color1 = list(cmd.get_color_tuple(color1))
    color2 = list(cmd.get_color_tuple(color2))

    def get_coord(v):
        if not isinstance(v, str):
            return v
        if v.startswith('['):
            return cmd.safe_list_eval(v)
        return cmd.get_atom_coords(v)

    xyz1 = get_coord(atom1)
    xyz2 = get_coord(atom2)
    normal = cpv.normalize(cpv.sub(xyz1, xyz2))

    if hlength < 0:
        hlength = radius * 3.0
    if hradius < 0:
        hradius = hlength * 0.6

    if gap:
        diff = cpv.scale(normal, gap)
        xyz1 = cpv.sub(xyz1, diff)
        xyz2 = cpv.add(xyz2, diff)

    xyz3 = cpv.add(cpv.scale(normal, hlength), xyz2)

    obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2 + \
          [cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 + \
          [1.0, 0.0]

    if not name:
        name = cmd.get_unused_name('arrow')

    cmd.load_cgo(obj, name)

cmd.extend('cgo_arrow', cgo_arrow)

def color_waters_yellow(): 
    yellow_waters = np.load("/Users/yyklab/Downloads/yellow_water_indices.npy")
    existing_selections = cmd.get_names("selections")

    #------COLOR YELLOW--------
    for water in yellow_waters: 
        sel_name = f"water_{water}"
        
        if sel_name in existing_selections:
            cmd.color("yellow", sel_name)
            cmd.set("label_color", "white", sel_name)
            

def main():

    cmd.bg_color("white")


    #----SELECTIONS-----
    cmd.load(protein_file)
    cmd.load(water_ccp4_map)
    cmd.select("chainG", "chain G")
    cmd.select("near_chainG", "(all within 5 of chainG) and not chainG")
    cmd.select("ASN_A_OD1_181", "chain A and resn ASN and resi 181")

    



    #----ISOLATION & GRAPHICS-----
    


    selections = []          
    for i in res: 
        parts = i.strip("'").split('_')
        
        selection = "chain A and resn " + parts[0] + " and resi " + parts[-1]
        print(selection)
        cmd.select(i, selection)   
        selections.append(selection)   

    cmd.hide("everything", "all")
    cmd.show("cartoon", "chainG")
    cmd.show("sticks", "near_chainG")

    selection = "near_chainG"

    for i in res: 
        selection += " and not " + i
    
    cmd.hide("sticks", selection)
    
    for i in res:
        if i.startswith("OEC"):
            cmd.show("spheres", i)
            cmd.set("sphere_scale", 0.25, i) 

            cmd.distance("OEC_sticks_"+ i, i, i, cutoff=2.6)
            cmd.label(i, "name")
            cmd.hide("labels", "OEC_sticks_" + i)      
            cmd.set("dash_gap", 0)                
            cmd.set("dash_radius", 0.15)          
            cmd.set("dash_round_ends", 1)     
     
                   

    # Draw all arrows
    arrows = [('resn OOO and chain G and resi 48',
  'resn ASP and chain A and name O and resi 61',
  'red1'),
 ('resn OOO and chain G and resi 21',
  'resn ASP and chain A and name OD2 and resi 61',
  'red2'),
 ('resn OOO and chain G and resi 62',
  'resn SER and chain A and name O and resi 169',
  'red1'),
 ('resn OEC and chain A and name W1 and resi 601',
  'resn SER and chain A and name O and resi 169',
  'red1'),
 ('resn OOO and chain G and resi 62',
  'resn GLY and chain A and name O and resi 171',
  'red2'),
 ('resn ASN and chain A and name OD1 and resi 181',
  'resn ASN and chain A and name ND2 and resi 181',
  'red1'),
 ('resn OOO and chain G and resi 52',
  'resn ASN and chain A and name O and resi 335',
  'red3'),
 ('resn OOO and chain G and resi 34',
  'resn PRO and chain A and name O and resi 340',
  'red1'),
 ('resn OOO and chain G and resi 26',
  'resn OEC and chain A and name O1 and resi 601',
  'red1'),
 ('resn OOO and chain G and resi 19',
  'resn OEC and chain A and name O4 and resi 601',
  'red1'),
 ('resn OOO and chain G and resi 21',
  'resn OEC and chain A and name W1 and resi 601',
  'red2'),
 ('resn OOO and chain G and resi 22',
  'resn OEC and chain A and name W2 and resi 601',
  'red2'),
 ('resn OOO and chain G and resi 28',
  'resn OEC and chain A and name W4 and resi 601',
  'red3'),
 ('resn ARG and chain A and name NH1 and resi 334',
  'resn GLU and chain D and name OE1 and resi 312',
  'red2'),
 ('resn OOO and chain G and resi 66',
  'resn GLU and chain D and name OE1 and resi 312',
  'red1'),
 ('resn ARG and chain A and name NH1 and resi 334',
  'resn GLU and chain D and name OE2 and resi 312',
  'red2'),
 ('resn OOO and chain G and resi 66',
  'resn GLU and chain D and name OE2 and resi 312',
  'red1'),
 ('resn ASN and chain A and name ND2 and resi 181',
  'resn OOO and chain G and resi 22',
  'red2'),
 ('resn OOO and chain G and resi 24',
  'resn OOO and chain G and resi 23',
  'red2'),
 ('resn OOO and chain G and resi 61',
  'resn OOO and chain G and resi 24',
  'red2'),
 ('resn OEC and chain A and name W2 and resi 601',
  'resn OOO and chain G and resi 24',
  'red2'),
 ('resn ASN and chain A and name ND2 and resi 181',
  'resn OOO and chain G and resi 61',
  'red1'),
 ('resn OEC and chain A and name W2 and resi 601',
  'resn OOO and chain G and resi 61',
  'red3'),
 ('resn OOO and chain G and resi 61',
  'resn OOO and chain G and resi 62',
  'red1'),
 ('resn GLU and chain A and name OE1 and resi 65',
  'resn OOO and chain G and resi 42',
  'red3'),
 ('resn GLU and chain A and name OE2 and resi 65',
  'resn OOO and chain G and resi 42',
  'red3')]

    for i, (start, end, color) in enumerate(arrows):
        cgo_arrow(start, end, color=color, radius=0.2, hlength=0.4, hradius=0.4, gap=0.6, name=f'arrow_{i+1}')
    #----PEAK PDB FILES-----
    cmd.load(H1_peaks)
    cmd.load(H2_peaks)
    cmd.load(OW_peaks)

    #----SHOW SPHERES-----
    spheres = ['show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`1551/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`1617/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`399/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`280/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`4900/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`480/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`897/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`28/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`135/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`2668/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`921/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`4035/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`1538/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`1538/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`2713/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`6345/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`2919/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`5358/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`462/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`3360/O',
 'show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`1958/O',
 'show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`1100/O',
 'show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`1329/O',
 'show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`1108/O',
 'show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`1411/O',
 'show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`1909/O',
 'show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`2692/O',
 'show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`702/O',
 'show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`1644/O',
 'show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`2355/O',
 'show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`2704/O',
 'show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`2691/O',
 'show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`2918/O',
 'show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`4974/O',
 'show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`3233/O',
 'show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`4016/O',
 'show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`1603/O',
 'show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`3412/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`462/O',
 'show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`702/O',
 'show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`1603/O',
 'show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`263/O',
 'show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`976/O',
 'show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`1466/O']
    for sphere in spheres:
        cmd.show("spheres", sphere[13:])

    #----HIDE EXTRA PEAKS-----
    cmd.hide("nonbonded", "peaks_low_restraints_justH1_1_rms_2.13")
    cmd.hide("nonbonded", "peaks_low_restraints_justH2_1_rms_2.11")
    cmd.hide("nonbonded", "peaks_low_restraints_justOW_1_rms_1.58")


    #----SELECT WATERS-----
    protein_stem = Path(protein_file).stem
    for resi in waters:
        individual_name = f"water_{resi}"
        selection_string = f"/{protein_stem}//G/{resi}/O"
        
        cmd.select(individual_name, selection_string)
        
        cmd.show("spheres", individual_name)
        cmd.set("sphere_scale", 0.25, individual_name) 
        cmd.color("blue", individual_name)
        
        cmd.label(individual_name, f'"W{resi}"')
        cmd.show("labels", individual_name)        
        cmd.set("label_size", 14, individual_name)  
        cmd.set("label_color", "white")
        cmd.set("label_font_id", 7)
    #----COLOR YELLOW-----
    color_waters_yellow()

    cmd.isomesh("meshmap", "fcalc_lower_restraints_last55ns_p212121_justwater_1", 1.0, selection = "0F.pdb_fittedto55all and chain G", carve = 2.0)
    cmd.color_deep("gray40", 'meshmap', 0)






main()