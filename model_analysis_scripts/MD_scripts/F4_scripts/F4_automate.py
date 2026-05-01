
import numpy as np

pickle_stem = "/Users/yyklab/Downloads"
file_stem = "/Users/yyklab/Desktop/Lab_Files/low_restraints"

#----INPUTS-----
protein_file = "/Users/yyklab/Desktop/Lab_Files/restart_09_09_25/base_files/0F.pdb_fittedto55all.pdb" # 0F.pdb
H1_peaks = "/Users/yyklab/Desktop/Lab_Files/low_restraints/maps_peaks/peaks_low_restraints_justH1_1_rms_2.13.pdb"
H2_peaks = "/Users/yyklab/Desktop/Lab_Files/low_restraints/maps_peaks/peaks_low_restraints_justH2_1_rms_2.11.pdb"
OW_peaks = "/Users/yyklab/Desktop/Lab_Files/low_restraints/maps_peaks/peaks_low_restraints_justOW_1_rms_1.58.pdb"
water_ccp4_map = "/Users/yyklab/Desktop/Lab_Files/low_restraints/base_files/FFT_6/fcalc_lower_restraints_last55ns_p212121_justwater_1.ccp4"

def O1(chain, finalgro, st_O1, ax1, ax2, donor, acceptor, strength, ticklabels1, ticklabels2, OW_peaks_involved, H1_peaks_involved, H2_peaks_involved, H1_partners_involved, H2_partners_involved, plot_dub, chan):
    import pandas as pd
    import math
    import matplotlib.pyplot as plt
    import matplotlib.ticker as plticker
    import matplotlib.patches as mpatches
    import matplotlib
    from IPython.display import clear_output
    import timeit
    from pathlib import Path
    if chain == 'G':
        offset = 0
    else:
        offset = 147
        
    haystack= list((finalgro['resseq'][offset:offset+147]).values)


    o1_b =[]
    for x in st_O1:
        if x in haystack:
            o1_b.append((haystack.index(x)+offset))
    z=0
    prev1=0
    prev2 = 0 

    for p in o1_b:
        # Skip if no h partners

        if (finalgro['waterID_h2_partner_fin'].iloc[p] ==0 and finalgro['waterID_h1_partner_fin'].iloc[p] ==0 ):
        
            continue
        #Skip if no oxygen
        if finalgro['OW_ID'].iloc[p]==0:
            continue
        OW_peaks_involved.append(finalgro['OW_ID'].iloc[p])
        
        H1_peaks_involved.append(finalgro['HW1_ID'].iloc[p])
        H2_peaks_involved.append(finalgro['HW2_ID'].iloc[p])  


        if (finalgro['waterID_h1_partner_fin'].iloc[p] == finalgro['waterID_h2_partner_fin'].iloc[p]) and finalgro['ow_h1_partner_dist_fin'].iloc[p] < finalgro['ow_h2_partner_dist'].iloc[p]:
            strength.append(finalgro['ow_h1_partner_dist_fin'].iloc[p])
        
            ax1.scatter(z,finalgro['ow_h1_partner_dist_fin'].iloc[p],s=100, color='blue',edgecolor='black',label ='O-H1 distance')
            ax1.scatter(z,finalgro['o_o_dist_h1_partner_fin'].iloc[p] ,s=100,color='orange',edgecolor='black', label ='O-O distance')
            ax2.plot(z, 1,color='white') 
            ticklabels1.append(finalgro['resseq'].iloc[p])
            ticklabels2.append(finalgro['waterID_h1_partner_fin'].iloc[p])
            
            prev1 = finalgro['waterID_h1_partner_fin'].iloc[p]
        elif (finalgro['waterID_h1_partner'].iloc[p] == finalgro['waterID_h2_partner_fin'].iloc[p]) and finalgro['ow_h2_partner_dist_fin'].iloc[p] < finalgro['ow_h1_partner_dist'].iloc[p]:
            strength.append(finalgro['ow_h2_partner_dist_fin'].iloc[p])
            ax1.scatter(z,finalgro['ow_h2_partner_dist_fin'].iloc[p],s=100, color='blue',edgecolor='black',label ='O-H1 distance')
            ax1.scatter(z,finalgro['o_o_dist_h2_partner_fin'].iloc[p] ,s=100,color='orange',edgecolor='black', label ='O-O distance')
            ax2.plot(z, 1,color='white') 
            ticklabels1.append(finalgro['resseq'].iloc[p])
            ticklabels2.append(finalgro['waterID_h2_partner_fin'].iloc[p])
            prev1 = finalgro['waterID_h2_partner_fin'].iloc[p]
        
        elif (finalgro['waterID_h2_partner_fin'].iloc[p]==0): 
            strength.append(finalgro['ow_h1_partner_dist_fin'].iloc[p])
            ax1.scatter(z,finalgro['ow_h1_partner_dist_fin'].iloc[p],s=100, color='blue',edgecolor='black',label ='O-H1 distance')
            ax1.scatter(z,finalgro['o_o_dist_h1_partner_fin'].iloc[p] ,s=100,color='orange',edgecolor='black', label ='O-O distance')
            ax2.plot(z, 1,color='white') 
            ticklabels1.append(finalgro['resseq'].iloc[p])
            ticklabels2.append(finalgro['waterID_h1_partner_fin'].iloc[p])
            prev1 = finalgro['waterID_h1_partner_fin'].iloc[p]
        elif (finalgro['waterID_h1_partner_fin'].iloc[p]==0): 
            strength.append(finalgro['ow_h2_partner_dist_fin'].iloc[p])
            ax1.scatter(z,finalgro['ow_h2_partner_dist_fin'].iloc[p],s=100, color='blue',edgecolor='black',label ='O-H1 distance')
            ax1.scatter(z,finalgro['o_o_dist_h2_partner_fin'].iloc[p] ,s=100,color='orange',edgecolor='black', label ='O-O distance')
            ax2.plot(z, 1,color='white') 
            ticklabels1.append(finalgro['resseq'].iloc[p])
            ticklabels2.append(finalgro['waterID_h2_partner_fin'].iloc[p])
            prev1 = finalgro['waterID_h2_partner_fin'].iloc[p]
        else:
            z, prev1 = plot_dub(z,p,chain)
        
        H2_partners_involved.append(finalgro['waterID_h2_partner_fin'].iloc[p])
        H1_partners_involved.append(finalgro['waterID_h1_partner_fin'].iloc[p])   
        
        z=z+1
        


    ax1.set_xticks(np.arange(0,len(ticklabels2),1))
    ax2.set_xticks(np.arange(0,len(ticklabels2),1))

    # NORMALIZATION FOR DIFF NUMPY VERSIONS
    def normalize_label(x):
        if isinstance(x, (list, tuple)):
            return '' if len(x) == 0 else str(x[0])  # or join with ',' if multiple
        return str(x)

    ticklabels2 = [normalize_label(x) for x in ticklabels2]
    print(ticklabels2)


    ax2.set_xticklabels(np.array(ticklabels2))
    ax1.set_xticklabels(np.array(ticklabels1))

    donor.append(ticklabels2)
    acceptor.append(ticklabels1)
    ax1.fill_between(np.arange(0,len(ticklabels2),1), 2.6, 3, color='lightpink',alpha=0.1, label ='weak')
    ax1.fill_between(np.arange(0,len(ticklabels2),1), 2.4, 2.6, color='salmon',alpha = 0.3, label='moderate')
    ax1.fill_between(np.arange(0,len(ticklabels2),1), 0, 2.4, color='red', alpha = 0.4, label='strong')

    a = mpatches.Patch(color='lightpink', alpha = 0.1,  label ='weak')
    b= mpatches.Patch(color='salmon', alpha = 0.3,  label='moderate')
    c= mpatches.Patch(color='red', alpha = 0.4,  label='strong')
    d= mpatches.Patch(color='blue',  label='o-h distance')
    e= mpatches.Patch(color='orange',  label='o-o distance')

    ax2.legend(handles=[a,b,c], bbox_to_anchor=(1.1, 0.8))

    colors = ["b", "orange"]
    texts = ["o-h distance", "o-o distance"]
    patches = [ plt.plot([],[], marker="o", ms=10, ls="", mec=None, color=colors[i], 
                label="{:s}".format(texts[i]) )[0]  for i in range(len(texts)) ]
    ax1.legend(handles=patches, bbox_to_anchor=(1, 1), 
            loc='upper left', ncol=2,  numpoints=1)
    ax1.set_ylim(1,3.5)
    ax1.set_ylim(1,3.5)
    loc = plticker.MultipleLocator(base=1.0) # this locator puts ticks at regular intervals
    plt.tight_layout()
    ax1.tick_params(labelrotation=90)
    ax2.tick_params(labelrotation=90)
    if chain =='G':
        label ='firstchain'
    else:
        label='secondchain'
    plt.savefig('Kanishk_%s_hbond_%s_twoH.png' %(chan, label),  bbox_inches='tight',dpi=100)

    # cell 26

    finalgro[
        (finalgro['chainid'] == 'G') &
        (finalgro['resseq'].isin(range(20, 31)))
    ]

    # cell 27

    arrows = []
    df = pd.DataFrame()
    df['donor'] = donor[0]
    df['acceptor'] = acceptor[0]
    df['strength'] = strength

    df.to_pickle(Path(file_stem) / ('%s_%s_donor_acceptor_strength_pickle' %(chan, label)))
    example=pd.read_pickle(Path(file_stem) / ('%s_%s_donor_acceptor_strength_pickle' %(chan, label)))


        
    for z in range(0,len(donor[0])):
        if strength[z] <= 2.4:
            color= 'red1' # red 1
        elif strength[z] > 2.4 and strength[z] <= 2.6 :
            color='red2' # red 2
        elif strength[z] > 2.6 and strength[z] <= 3 :
            color='red3' # red 3
        try:
            float(example.loc[z,'donor']) #First assume you are only dealing with chain G waters
            don = 'resn OOO and chain %s and resi %d' %(chain,float(example.loc[z,'donor']))
        except ValueError:
            new_name = example['donor'][z].split('_')
            don = 'resn %s and chain %s and name %s and resi %s' %(new_name[0], new_name[1], new_name[2], new_name[3])
        
                
        try:
            float(example.loc[z,'acceptor'])
            acc = 'resn OOO and chain %s and resi %d' %(chain,float(example.loc[z,'acceptor']))
        except ValueError:
            new_name = example['acceptor'][z].split('_')
            acc = 'resn %s and chain %s and name %s and resi %s' %(new_name[0], new_name[1], new_name[2], new_name[3])

        don = don.replace("'", '')  
        acc = acc.replace("'", '') 

        arrows.append((don, acc, color))
        #arrows.append('(%s, %s,  color=%s, radius= 0.2, hlength=0.4, hradius=0.4, gap = 0.6)' %(don, acc,color))

    print(arrows)
    # cell 28

    a = 0
    spheres = []
    spheres.append
    for z in OW_peaks_involved:
        if z ==0:
            continue
        else:
            spheres.append('show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`%d/O' %(z))
            

    for z in H1_peaks_involved:
        if z ==0:
            continue
        else:
            spheres.append('show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`%d/O' %(z))
                    
                    
    for z in H2_peaks_involved:
        if z ==0:
            continue
        else:
            spheres.append('show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`%d/O' %(z))

        


    for z in H1_partners_involved:
        if z ==0:
            continue
        #FIX
        if isinstance(z, list):
            if len(z) == 0:
                continue
            z = z[0]
            
        else:
            try:
                float(z)
                
                ow = finalgro[finalgro['resseq'] == float(z)]
                h1 = finalgro[finalgro['resseq'] == float(z)]
                h2 = finalgro[finalgro['resseq'] == float(z)]
                
                if len(ow) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`%d/O' %(ow.iloc[0]['OW_ID']))
                
                if len(h1) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`%d/O' %(h1.iloc[0]['HW1_ID']))
                
                if len(h2) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`%d/O' %(h2.iloc[0]['HW2_ID']))
                
            except (ValueError, TypeError): 
    
                ow = finalgro[finalgro['resseq'] == z.replace("'", '')]
                h1 = finalgro[finalgro['resseq'] == z.replace("'", '')]
                h2 = finalgro[finalgro['resseq'] == z.replace("'", '')]

                if len(ow) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`%d/O' %(ow.iloc[0]['OW_ID']))
                
                if len(h1) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`%d/O' %(h1.iloc[0]['HW1_ID']))
                
                if len(h2) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`%d/O' %(h2.iloc[0]['HW2_ID']))
                else:
                    continue
                    
    for z in H2_partners_involved:
        if z ==0:
            continue
        #FIX
        if isinstance(z, list):
            if len(z) == 0:
                continue
            z = z[0]
        #FIX
        else:
            try:
                float(z)
                
                ow = finalgro[finalgro['resseq'] == float(z)]
                h1 = finalgro[finalgro['resseq'] == float(z)]
                h2 = finalgro[finalgro['resseq'] == float(z)]
                
                if len(ow) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`%d/O' %(ow.iloc[0]['OW_ID']))
                
                if len(h1) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`%d/O' %(h1.iloc[0]['HW1_ID']))
                
                if len(h2) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`%d/O' %(h2.iloc[0]['HW2_ID']))
                
            except ValueError: 

    
                ow = finalgro[finalgro['resseq'] == z.replace("'", '') ]
                h1 = finalgro[finalgro['resseq'] == z.replace("'", '')]
                h2 = finalgro[finalgro['resseq'] == z.replace("'", '')]

                if len(ow) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`%d/O' %(ow.iloc[0]['OW_ID']))
                
                if len(h1) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`%d/O' %(h1.iloc[0]['HW1_ID']))
                
                if len(h2) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`%d/O' %(h2.iloc[0]['HW2_ID']))
                else:
                    continue
    np.save(Path(pickle_stem) / 'ticklabels1_O1',ticklabels1)
    np.save(Path(pickle_stem) / 'ticklabels2_O1',ticklabels2)
    print(spheres)
    np.save(Path(pickle_stem) / 'spheres_O1', spheres)
    np.save(Path(pickle_stem) / 'arrows_O1', arrows)

    # cell 29

def Cl(chain, finalgro, st_Cl, ax1, ax2, donor, acceptor, strength, ticklabels1, ticklabels2, OW_peaks_involved, H1_peaks_involved, H2_peaks_involved, H1_partners_involved, H2_partners_involved, plot_dub, chan):
    import pandas as pd
    import math
    import matplotlib.pyplot as plt
    import matplotlib.ticker as plticker
    import matplotlib.patches as mpatches
    import matplotlib
    from IPython.display import clear_output
    import timeit
    from pathlib import Path

    if chain == 'G':
        offset = 0
    else:
        offset = 147
        
    haystack= list((finalgro['resseq'][offset:offset+147]).values)


    o1_b =[]
    for x in st_Cl:
        if x in haystack:
            o1_b.append((haystack.index(x)+offset))
    z=0
    prev1=0
    prev2 = 0 

    for p in o1_b:
        # Skip if no h partners

        if (finalgro['waterID_h2_partner_fin'].iloc[p] ==0 and finalgro['waterID_h1_partner_fin'].iloc[p] ==0 ):
        
            continue
        #Skip if no oxygen
        if finalgro['OW_ID'].iloc[p]==0:
            continue
        OW_peaks_involved.append(finalgro['OW_ID'].iloc[p])
        
        H1_peaks_involved.append(finalgro['HW1_ID'].iloc[p])
        H2_peaks_involved.append(finalgro['HW2_ID'].iloc[p])  


        if (finalgro['waterID_h1_partner_fin'].iloc[p] == finalgro['waterID_h2_partner_fin'].iloc[p]) and finalgro['ow_h1_partner_dist_fin'].iloc[p] < finalgro['ow_h2_partner_dist'].iloc[p]:
            strength.append(finalgro['ow_h1_partner_dist_fin'].iloc[p])
        
            ax1.scatter(z,finalgro['ow_h1_partner_dist_fin'].iloc[p],s=100, color='blue',edgecolor='black',label ='O-H1 distance')
            ax1.scatter(z,finalgro['o_o_dist_h1_partner_fin'].iloc[p] ,s=100,color='orange',edgecolor='black', label ='O-O distance')
            ax2.plot(z, 1,color='white') 
            ticklabels1.append(finalgro['resseq'].iloc[p])
            ticklabels2.append(finalgro['waterID_h1_partner_fin'].iloc[p])
            
            prev1 = finalgro['waterID_h1_partner_fin'].iloc[p]
        elif (finalgro['waterID_h1_partner'].iloc[p] == finalgro['waterID_h2_partner_fin'].iloc[p]) and finalgro['ow_h2_partner_dist_fin'].iloc[p] < finalgro['ow_h1_partner_dist'].iloc[p]:
            strength.append(finalgro['ow_h2_partner_dist_fin'].iloc[p])
            ax1.scatter(z,finalgro['ow_h2_partner_dist_fin'].iloc[p],s=100, color='blue',edgecolor='black',label ='O-H1 distance')
            ax1.scatter(z,finalgro['o_o_dist_h2_partner_fin'].iloc[p] ,s=100,color='orange',edgecolor='black', label ='O-O distance')
            ax2.plot(z, 1,color='white') 
            ticklabels1.append(finalgro['resseq'].iloc[p])
            ticklabels2.append(finalgro['waterID_h2_partner_fin'].iloc[p])
            prev1 = finalgro['waterID_h2_partner_fin'].iloc[p]
        
        elif (finalgro['waterID_h2_partner_fin'].iloc[p]==0): 
            strength.append(finalgro['ow_h1_partner_dist_fin'].iloc[p])
            ax1.scatter(z,finalgro['ow_h1_partner_dist_fin'].iloc[p],s=100, color='blue',edgecolor='black',label ='O-H1 distance')
            ax1.scatter(z,finalgro['o_o_dist_h1_partner_fin'].iloc[p] ,s=100,color='orange',edgecolor='black', label ='O-O distance')
            ax2.plot(z, 1,color='white') 
            ticklabels1.append(finalgro['resseq'].iloc[p])
            ticklabels2.append(finalgro['waterID_h1_partner_fin'].iloc[p])
            prev1 = finalgro['waterID_h1_partner_fin'].iloc[p]
        elif (finalgro['waterID_h1_partner_fin'].iloc[p]==0): 
            strength.append(finalgro['ow_h2_partner_dist_fin'].iloc[p])
            ax1.scatter(z,finalgro['ow_h2_partner_dist_fin'].iloc[p],s=100, color='blue',edgecolor='black',label ='O-H1 distance')
            ax1.scatter(z,finalgro['o_o_dist_h2_partner_fin'].iloc[p] ,s=100,color='orange',edgecolor='black', label ='O-O distance')
            ax2.plot(z, 1,color='white') 
            ticklabels1.append(finalgro['resseq'].iloc[p])
            ticklabels2.append(finalgro['waterID_h2_partner_fin'].iloc[p])
            prev1 = finalgro['waterID_h2_partner_fin'].iloc[p]
        else:
            z,prev1 = plot_dub(z,p,chain)
        
        H2_partners_involved.append(finalgro['waterID_h2_partner_fin'].iloc[p])
        H1_partners_involved.append(finalgro['waterID_h1_partner_fin'].iloc[p])   
        
        z=z+1
        


    ax1.set_xticks(np.arange(0,len(ticklabels2),1))
    ax2.set_xticks(np.arange(0,len(ticklabels2),1))

    # NORMALIZATION FOR DIFF NUMPY VERSIONS
    def normalize_label(x):
        if isinstance(x, (list, tuple)):
            return '' if len(x) == 0 else str(x[0])  # or join with ',' if multiple
        return str(x)

    ticklabels2 = [normalize_label(x) for x in ticklabels2]
    print(ticklabels2)


    ax2.set_xticklabels(np.array(ticklabels2))
    ax1.set_xticklabels(np.array(ticklabels1))

    donor.append(ticklabels2)
    acceptor.append(ticklabels1)
    ax1.fill_between(np.arange(0,len(ticklabels2),1), 2.6, 3, color='lightpink',alpha=0.1, label ='weak')
    ax1.fill_between(np.arange(0,len(ticklabels2),1), 2.4, 2.6, color='salmon',alpha = 0.3, label='moderate')
    ax1.fill_between(np.arange(0,len(ticklabels2),1), 0, 2.4, color='red', alpha = 0.4, label='strong')

    a = mpatches.Patch(color='lightpink', alpha = 0.1,  label ='weak')
    b= mpatches.Patch(color='salmon', alpha = 0.3,  label='moderate')
    c= mpatches.Patch(color='red', alpha = 0.4,  label='strong')
    d= mpatches.Patch(color='blue',  label='o-h distance')
    e= mpatches.Patch(color='orange',  label='o-o distance')

    ax2.legend(handles=[a,b,c], bbox_to_anchor=(1.1, 0.8))

    colors = ["b", "orange"]
    texts = ["o-h distance", "o-o distance"]
    patches = [ plt.plot([],[], marker="o", ms=10, ls="", mec=None, color=colors[i], 
                label="{:s}".format(texts[i]) )[0]  for i in range(len(texts)) ]
    ax1.legend(handles=patches, bbox_to_anchor=(1, 1), 
            loc='upper left', ncol=2,  numpoints=1)
    ax1.set_ylim(1,3.5)
    ax1.set_ylim(1,3.5)
    loc = plticker.MultipleLocator(base=1.0) # this locator puts ticks at regular intervals
    plt.tight_layout()
    ax1.tick_params(labelrotation=90)
    ax2.tick_params(labelrotation=90)
    if chain =='G':
        label ='firstchain'
    else:
        label='secondchain'
    plt.savefig('Kanishk_%s_hbond_%s_twoH.png' %(chan, label),  bbox_inches='tight',dpi=100)

    # cell 26

    finalgro[
        (finalgro['chainid'] == 'G') &
        (finalgro['resseq'].isin(range(20, 31)))
    ]

    # cell 27

    arrows = []
    df = pd.DataFrame()
    df['donor'] = donor[0]
    df['acceptor'] = acceptor[0]
    df['strength'] = strength

    df.to_pickle(Path(file_stem) / ('%s_%s_donor_acceptor_strength_pickle' %(chan, label)))
    example=pd.read_pickle(Path(file_stem) / ('%s_%s_donor_acceptor_strength_pickle' %(chan, label)))


        
    for z in range(0,len(donor[0])):
        if strength[z] <= 2.4:
            color= 'red1' # red 1
        elif strength[z] > 2.4 and strength[z] <= 2.6 :
            color='red2' # red 2
        elif strength[z] > 2.6 and strength[z] <= 3 :
            color='red3' # red 3
        try:
            float(example.loc[z,'donor']) #First assume you are only dealing with chain G waters
            don = 'resn OOO and chain %s and resi %d' %(chain,float(example.loc[z,'donor']))
        except ValueError:
            new_name = example['donor'][z].split('_')
            don = 'resn %s and chain %s and name %s and resi %s' %(new_name[0], new_name[1], new_name[2], new_name[3])
        
                
        try:
            float(example.loc[z,'acceptor'])
            acc = 'resn OOO and chain %s and resi %d' %(chain,float(example.loc[z,'acceptor']))
        except ValueError:
            new_name = example['acceptor'][z].split('_')
            acc = 'resn %s and chain %s and name %s and resi %s' %(new_name[0], new_name[1], new_name[2], new_name[3])

        don = don.replace("'", '')  
        acc = acc.replace("'", '') 

        arrows.append((don, acc, color))
        #arrows.append('(%s, %s,  color=%s, radius= 0.2, hlength=0.4, hradius=0.4, gap = 0.6)' %(don, acc,color))

    print(arrows)
    # cell 28

    a = 0
    spheres = []
    spheres.append
    for z in OW_peaks_involved:
        if z ==0:
            continue
        else:
            spheres.append('show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`%d/O' %(z))
            

    for z in H1_peaks_involved:
        if z ==0:
            continue
        else:
            spheres.append('show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`%d/O' %(z))
                    
                    
    for z in H2_peaks_involved:
        if z ==0:
            continue
        else:
            spheres.append('show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`%d/O' %(z))

        


    for z in H1_partners_involved:
        if z ==0:
            continue
        #FIX
        if isinstance(z, list):
            if len(z) == 0:
                continue
            z = z[0]
            
        else:
            try:
                float(z)
                
                ow = finalgro[finalgro['resseq'] == float(z)]
                h1 = finalgro[finalgro['resseq'] == float(z)]
                h2 = finalgro[finalgro['resseq'] == float(z)]
                
                if len(ow) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`%d/O' %(ow.iloc[0]['OW_ID']))
                
                if len(h1) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`%d/O' %(h1.iloc[0]['HW1_ID']))
                
                if len(h2) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`%d/O' %(h2.iloc[0]['HW2_ID']))
                
            except (ValueError, TypeError): 
    
                ow = finalgro[finalgro['resseq'] == z.replace("'", '')]
                h1 = finalgro[finalgro['resseq'] == z.replace("'", '')]
                h2 = finalgro[finalgro['resseq'] == z.replace("'", '')]

                if len(ow) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`%d/O' %(ow.iloc[0]['OW_ID']))
                
                if len(h1) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`%d/O' %(h1.iloc[0]['HW1_ID']))
                
                if len(h2) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`%d/O' %(h2.iloc[0]['HW2_ID']))
                else:
                    continue
                    
    for z in H2_partners_involved:
        if z ==0:
            continue
        #FIX
        if isinstance(z, list):
            if len(z) == 0:
                continue
            z = z[0]
        #FIX
        else:
            try:
                float(z)
                
                ow = finalgro[finalgro['resseq'] == float(z)]
                h1 = finalgro[finalgro['resseq'] == float(z)]
                h2 = finalgro[finalgro['resseq'] == float(z)]
                
                if len(ow) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`%d/O' %(ow.iloc[0]['OW_ID']))
                
                if len(h1) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`%d/O' %(h1.iloc[0]['HW1_ID']))
                
                if len(h2) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`%d/O' %(h2.iloc[0]['HW2_ID']))
                
            except ValueError: 

    
                ow = finalgro[finalgro['resseq'] == z.replace("'", '') ]
                h1 = finalgro[finalgro['resseq'] == z.replace("'", '')]
                h2 = finalgro[finalgro['resseq'] == z.replace("'", '')]

                if len(ow) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`%d/O' %(ow.iloc[0]['OW_ID']))
                
                if len(h1) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`%d/O' %(h1.iloc[0]['HW1_ID']))
                
                if len(h2) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`%d/O' %(h2.iloc[0]['HW2_ID']))
                else:
                    continue
    np.save(Path(pickle_stem) / 'ticklabels1_Cl',ticklabels1)
    np.save(Path(pickle_stem) / 'ticklabels2_Cl',ticklabels2)
    print(spheres)
    np.save(Path(pickle_stem) / 'spheres_Cl', spheres)
    np.save(Path(pickle_stem) / 'arrows_Cl', arrows)


def O4(chain, finalgro, st_O4, ax1, ax2, donor, acceptor, strength, ticklabels1, ticklabels2, OW_peaks_involved, H1_peaks_involved, H2_peaks_involved, H1_partners_involved, H2_partners_involved, plot_dub, chan):
    import pandas as pd
    import math
    import matplotlib.pyplot as plt
    import matplotlib.ticker as plticker
    import matplotlib.patches as mpatches
    import matplotlib
    from IPython.display import clear_output
    import timeit
    from pathlib import Path

    if chain == 'G':
        offset = 0
    else:
        offset = 147
        
    haystack= list((finalgro['resseq'][offset:offset+147]).values)


    o1_b =[]
    for x in st_O4:
        if x in haystack:
            o1_b.append((haystack.index(x)+offset))
    z=0
    prev1=0
    prev2 = 0 

    for p in o1_b:
        # Skip if no h partners

        if (finalgro['waterID_h2_partner_fin'].iloc[p] ==0 and finalgro['waterID_h1_partner_fin'].iloc[p] ==0 ):
        
            continue
        #Skip if no oxygen
        if finalgro['OW_ID'].iloc[p]==0:
            continue
        OW_peaks_involved.append(finalgro['OW_ID'].iloc[p])
        
        H1_peaks_involved.append(finalgro['HW1_ID'].iloc[p])
        H2_peaks_involved.append(finalgro['HW2_ID'].iloc[p])  


        if (finalgro['waterID_h1_partner_fin'].iloc[p] == finalgro['waterID_h2_partner_fin'].iloc[p]) and finalgro['ow_h1_partner_dist_fin'].iloc[p] < finalgro['ow_h2_partner_dist'].iloc[p]:
            strength.append(finalgro['ow_h1_partner_dist_fin'].iloc[p])
        
            ax1.scatter(z,finalgro['ow_h1_partner_dist_fin'].iloc[p],s=100, color='blue',edgecolor='black',label ='O-H1 distance')
            ax1.scatter(z,finalgro['o_o_dist_h1_partner_fin'].iloc[p] ,s=100,color='orange',edgecolor='black', label ='O-O distance')
            ax2.plot(z, 1,color='white') 
            ticklabels1.append(finalgro['resseq'].iloc[p])
            ticklabels2.append(finalgro['waterID_h1_partner_fin'].iloc[p])
            
            prev1 = finalgro['waterID_h1_partner_fin'].iloc[p]
        elif (finalgro['waterID_h1_partner'].iloc[p] == finalgro['waterID_h2_partner_fin'].iloc[p]) and finalgro['ow_h2_partner_dist_fin'].iloc[p] < finalgro['ow_h1_partner_dist'].iloc[p]:
            strength.append(finalgro['ow_h2_partner_dist_fin'].iloc[p])
            ax1.scatter(z,finalgro['ow_h2_partner_dist_fin'].iloc[p],s=100, color='blue',edgecolor='black',label ='O-H1 distance')
            ax1.scatter(z,finalgro['o_o_dist_h2_partner_fin'].iloc[p] ,s=100,color='orange',edgecolor='black', label ='O-O distance')
            ax2.plot(z, 1,color='white') 
            ticklabels1.append(finalgro['resseq'].iloc[p])
            ticklabels2.append(finalgro['waterID_h2_partner_fin'].iloc[p])
            prev1 = finalgro['waterID_h2_partner_fin'].iloc[p]
        
        elif (finalgro['waterID_h2_partner_fin'].iloc[p]==0): 
            strength.append(finalgro['ow_h1_partner_dist_fin'].iloc[p])
            ax1.scatter(z,finalgro['ow_h1_partner_dist_fin'].iloc[p],s=100, color='blue',edgecolor='black',label ='O-H1 distance')
            ax1.scatter(z,finalgro['o_o_dist_h1_partner_fin'].iloc[p] ,s=100,color='orange',edgecolor='black', label ='O-O distance')
            ax2.plot(z, 1,color='white') 
            ticklabels1.append(finalgro['resseq'].iloc[p])
            ticklabels2.append(finalgro['waterID_h1_partner_fin'].iloc[p])
            prev1 = finalgro['waterID_h1_partner_fin'].iloc[p]
        elif (finalgro['waterID_h1_partner_fin'].iloc[p]==0): 
            strength.append(finalgro['ow_h2_partner_dist_fin'].iloc[p])
            ax1.scatter(z,finalgro['ow_h2_partner_dist_fin'].iloc[p],s=100, color='blue',edgecolor='black',label ='O-H1 distance')
            ax1.scatter(z,finalgro['o_o_dist_h2_partner_fin'].iloc[p] ,s=100,color='orange',edgecolor='black', label ='O-O distance')
            ax2.plot(z, 1,color='white') 
            ticklabels1.append(finalgro['resseq'].iloc[p])
            ticklabels2.append(finalgro['waterID_h2_partner_fin'].iloc[p])
            prev1 = finalgro['waterID_h2_partner_fin'].iloc[p]
        else:
            z,prev1 = plot_dub(z,p,chain)
        
        H2_partners_involved.append(finalgro['waterID_h2_partner_fin'].iloc[p])
        H1_partners_involved.append(finalgro['waterID_h1_partner_fin'].iloc[p])   
        
        z=z+1
        


    ax1.set_xticks(np.arange(0,len(ticklabels2),1))
    ax2.set_xticks(np.arange(0,len(ticklabels2),1))

    # NORMALIZATION FOR DIFF NUMPY VERSIONS
    def normalize_label(x):
        if isinstance(x, (list, tuple)):
            return '' if len(x) == 0 else str(x[0])  # or join with ',' if multiple
        return str(x)

    ticklabels2 = [normalize_label(x) for x in ticklabels2]
    print(ticklabels2)


    ax2.set_xticklabels(np.array(ticklabels2))
    ax1.set_xticklabels(np.array(ticklabels1))

    donor.append(ticklabels2)
    acceptor.append(ticklabels1)
    ax1.fill_between(np.arange(0,len(ticklabels2),1), 2.6, 3, color='lightpink',alpha=0.1, label ='weak')
    ax1.fill_between(np.arange(0,len(ticklabels2),1), 2.4, 2.6, color='salmon',alpha = 0.3, label='moderate')
    ax1.fill_between(np.arange(0,len(ticklabels2),1), 0, 2.4, color='red', alpha = 0.4, label='strong')

    a = mpatches.Patch(color='lightpink', alpha = 0.1,  label ='weak')
    b= mpatches.Patch(color='salmon', alpha = 0.3,  label='moderate')
    c= mpatches.Patch(color='red', alpha = 0.4,  label='strong')
    d= mpatches.Patch(color='blue',  label='o-h distance')
    e= mpatches.Patch(color='orange',  label='o-o distance')

    ax2.legend(handles=[a,b,c], bbox_to_anchor=(1.1, 0.8))

    colors = ["b", "orange"]
    texts = ["o-h distance", "o-o distance"]
    patches = [ plt.plot([],[], marker="o", ms=10, ls="", mec=None, color=colors[i], 
                label="{:s}".format(texts[i]) )[0]  for i in range(len(texts)) ]
    ax1.legend(handles=patches, bbox_to_anchor=(1, 1), 
            loc='upper left', ncol=2,  numpoints=1)
    ax1.set_ylim(1,3.5)
    ax1.set_ylim(1,3.5)
    loc = plticker.MultipleLocator(base=1.0) # this locator puts ticks at regular intervals
    plt.tight_layout()
    ax1.tick_params(labelrotation=90)
    ax2.tick_params(labelrotation=90)
    if chain =='G':
        label ='firstchain'
    else:
        label='secondchain'
    plt.savefig('Kanishk_%s_hbond_%s_twoH.png' %(chan, label),  bbox_inches='tight',dpi=100)

    # cell 26

    finalgro[
        (finalgro['chainid'] == 'G') &
        (finalgro['resseq'].isin(range(20, 31)))
    ]

    # cell 27

    arrows = []
    df = pd.DataFrame()
    df['donor'] = donor[0]
    df['acceptor'] = acceptor[0]
    df['strength'] = strength

    df.to_pickle(Path(file_stem) / ('%s_%s_donor_acceptor_strength_pickle' %(chan, label)))
    example=pd.read_pickle(Path(file_stem) / ('%s_%s_donor_acceptor_strength_pickle' %(chan, label)))


        
    for z in range(0,len(donor[0])):
        if strength[z] <= 2.4:
            color= 'red1' # red 1
        elif strength[z] > 2.4 and strength[z] <= 2.6 :
            color='red2' # red 2
        elif strength[z] > 2.6 and strength[z] <= 3 :
            color='red3' # red 3
        try:
            float(example.loc[z,'donor']) #First assume you are only dealing with chain G waters
            don = 'resn OOO and chain %s and resi %d' %(chain,float(example.loc[z,'donor']))
        except ValueError:
            new_name = example['donor'][z].split('_')
            don = 'resn %s and chain %s and name %s and resi %s' %(new_name[0], new_name[1], new_name[2], new_name[3])
        
                
        try:
            float(example.loc[z,'acceptor'])
            acc = 'resn OOO and chain %s and resi %d' %(chain,float(example.loc[z,'acceptor']))
        except ValueError:
            new_name = example['acceptor'][z].split('_')
            acc = 'resn %s and chain %s and name %s and resi %s' %(new_name[0], new_name[1], new_name[2], new_name[3])

        don = don.replace("'", '')  
        acc = acc.replace("'", '') 

        arrows.append((don, acc, color))
        #arrows.append('(%s, %s,  color=%s, radius= 0.2, hlength=0.4, hradius=0.4, gap = 0.6)' %(don, acc,color))

    print(arrows)
    # cell 28

    a = 0
    spheres = []
    spheres.append
    for z in OW_peaks_involved:
        if z ==0:
            continue
        else:
            spheres.append('show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`%d/O' %(z))
            

    for z in H1_peaks_involved:
        if z ==0:
            continue
        else:
            spheres.append('show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`%d/O' %(z))
                    
                    
    for z in H2_peaks_involved:
        if z ==0:
            continue
        else:
            spheres.append('show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`%d/O' %(z))

        


    for z in H1_partners_involved:
        if z ==0:
            continue
        #FIX
        if isinstance(z, list):
            if len(z) == 0:
                continue
            z = z[0]
            
        else:
            try:
                float(z)
                
                ow = finalgro[finalgro['resseq'] == float(z)]
                h1 = finalgro[finalgro['resseq'] == float(z)]
                h2 = finalgro[finalgro['resseq'] == float(z)]
                
                if len(ow) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`%d/O' %(ow.iloc[0]['OW_ID']))
                
                if len(h1) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`%d/O' %(h1.iloc[0]['HW1_ID']))
                
                if len(h2) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`%d/O' %(h2.iloc[0]['HW2_ID']))
                
            except (ValueError, TypeError): 
    
                ow = finalgro[finalgro['resseq'] == z.replace("'", '')]
                h1 = finalgro[finalgro['resseq'] == z.replace("'", '')]
                h2 = finalgro[finalgro['resseq'] == z.replace("'", '')]

                if len(ow) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`%d/O' %(ow.iloc[0]['OW_ID']))
                
                if len(h1) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`%d/O' %(h1.iloc[0]['HW1_ID']))
                
                if len(h2) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`%d/O' %(h2.iloc[0]['HW2_ID']))
                else:
                    continue
                    
    for z in H2_partners_involved:
        if z ==0:
            continue
        #FIX
        if isinstance(z, list):
            if len(z) == 0:
                continue
            z = z[0]
        #FIX
        else:
            try:
                float(z)
                
                ow = finalgro[finalgro['resseq'] == float(z)]
                h1 = finalgro[finalgro['resseq'] == float(z)]
                h2 = finalgro[finalgro['resseq'] == float(z)]
                
                if len(ow) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`%d/O' %(ow.iloc[0]['OW_ID']))
                
                if len(h1) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`%d/O' %(h1.iloc[0]['HW1_ID']))
                
                if len(h2) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`%d/O' %(h2.iloc[0]['HW2_ID']))
                
            except ValueError: 

    
                ow = finalgro[finalgro['resseq'] == z.replace("'", '') ]
                h1 = finalgro[finalgro['resseq'] == z.replace("'", '')]
                h2 = finalgro[finalgro['resseq'] == z.replace("'", '')]

                if len(ow) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justOW_1_rms_1.58//X/HOH`%d/O' %(ow.iloc[0]['OW_ID']))
                
                if len(h1) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH1_1_rms_2.13//X/HOH`%d/O' %(h1.iloc[0]['HW1_ID']))
                
                if len(h2) > 0:
                    spheres.append('show sphere, /peaks_low_restraints_justH2_1_rms_2.11//X/HOH`%d/O' %(h2.iloc[0]['HW2_ID']))
                else:
                    continue
    np.save(Path(pickle_stem) / 'ticklabels1_O4',ticklabels1)
    np.save(Path(pickle_stem) / 'ticklabels2_O4',ticklabels2)
    print(spheres)
    np.save(Path(pickle_stem) / 'spheres_O4', spheres)
    np.save(Path(pickle_stem) / 'arrows_O4', arrows)

def MADI_analysis():
    import pandas as pd
    import math
    import matplotlib.pyplot as plt
    import matplotlib.ticker as plticker
    import matplotlib.patches as mpatches
    import matplotlib
    from IPython.display import clear_output
    import timeit
    from pathlib import Path
    F_pdb_path =  Path(file_stem) / "base_files/0F.pdb_fittedto55all.pdb"
    OW_pdb_path = Path(file_stem) / "maps_peaks/peaks_low_restraints_justOW_1_rms_1.58.pdb"
    H1_pdb_path = Path(file_stem) / "maps_peaks/peaks_low_restraints_justH1_1_rms_2.13.pdb"
    H2_pdb_path = Path(file_stem) / "maps_peaks/peaks_low_restraints_justH2_1_rms_2.11.pdb"

    def flatten_list(_2d_list):
        flat_list = []
        # Iterate through the outer list
        for element in _2d_list:
            if type(element) is list:
                # If the element is of type list, iterate through the sublist
                for item in element:
                    flat_list.append(item)
            else:
                flat_list.append(element)
        return flat_list

    matplotlib.rcParams['axes.labelsize'] = 22
    matplotlib.rcParams['xtick.labelsize'] = 20
    matplotlib.rcParams['ytick.labelsize'] = 22
    matplotlib.rcParams['legend.fontsize'] = 22

    def distance(p0, p1):
        return math.sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2 +(p0[2] - p1[2])**2  )
    #cell 1

    hbond_res = pd.read_csv(Path(pickle_stem) / "hbond_res_updated.csv")
    len(hbond_res)

    # cell 2

    res_nums =[]
    O_sele=[]
    res_names=[]

    for i, rew in hbond_res.iterrows():
        O_sele.append(hbond_res['O selection'].iloc[i].split('/'))
        res_nums.append(hbond_res['Number'].iloc[i])
        res_names.append(hbond_res['RESN'].iloc[i])
        
    print("(reference['resname']=='OOO') | (reference['name']=='W1')| (reference['name']=='W2')|(reference['name']=='W3')| (reference['name']=='W4')|(reference['name']=='O1') & (reference['resname']=='OEC')|(reference['name']=='O2') & (reference['resname']=='OEC')|(reference['name']=='O3') & (reference['resname']=='OEC')|(reference['name']=='O4') & (reference['resname']=='OEC')|(reference['name']=='O5') & (reference['resname']=='OEC')|")

    a = 0
    ind = 0 
    for x in (O_sele):
        for z in x:
            a=a+1
            print('(reference["name"]== "%s") & (reference["resname"]== "%s") & (reference["resseq"]== %s)|' %(z, res_names[ind], res_nums[ind]), end='')
        ind +=1

    #cell 3

    colspecs_pdb = [(0, 6), (6, 11), (12, 16), (16, 17), (17, 20), (21, 22), (22, 26),
                (26, 27), (30, 38), (38, 46), (46, 54), (54, 60), (60, 66), (76, 78),
                (78, 80)]

    names_pdb = ['ATOM', 'serial', 'name', 'altloc', 'resname', 'chainid', 'resseq',
            'icode', 'x', 'y', 'z', 'occupancy', 'tempfactor', 'element', 'charge']

    reference = pd.read_fwf(F_pdb_path, names=names_pdb, colspecs=colspecs_pdb, skiprows=2)

    #cell 4

    filtered_ref=reference[(reference['resname']=='OOO') | (reference['name']=='W1')| (reference['name']=='W2')|(reference['name']=='W3')| (reference['name']=='W4')|(reference['name']=='O1') & (reference['resname']=='OEC')|(reference['name']=='O2') & (reference['resname']=='OEC')|(reference['name']=='O3') & (reference['resname']=='OEC')|(reference['name']=='O4') & (reference['resname']=='OEC')|(reference['name']=='O5') & (reference['resname']=='OEC')|
    (reference["name"]== "OG1") & (reference["resname"]== "THR") & (reference["resseq"]== 335)|(reference["name"]== "NH1") & (reference["resname"]== "ARG") & (reference["resseq"]== 334)|(reference["name"]== "NH2") & (reference["resname"]== "ARG") & (reference["resseq"]== 334)|(reference["name"]== "OE1") & (reference["resname"]== "GLU") & (reference["resseq"]== 65)|(reference["name"]== "OE2") & (reference["resname"]== "GLU") & (reference["resseq"]== 65)|(reference["name"]== "O") & (reference["resname"]== "PRO") & (reference["resseq"]== 340)|(reference["name"]== "OD1") & (reference["resname"]== "ASP") & (reference["resseq"]== 61)|(reference["name"]== "OD2") & (reference["resname"]== "ASP") & (reference["resseq"]== 61)|(reference["name"]== "O") & (reference["resname"]== "ASP") & (reference["resseq"]== 61)|(reference["name"]== "ND2") & (reference["resname"]== "ASN") & (reference["resseq"]== 181)|(reference["name"]== "OD1") & (reference["resname"]== "ASN") & (reference["resseq"]== 181)|(reference["name"]== "O") & (reference["resname"]== "ASN") & (reference["resseq"]== 181)|(reference["name"]== "O") & (reference["resname"]== "GLY") & (reference["resseq"]== 171)|(reference["name"]== "O") & (reference["resname"]== "SER") & (reference["resseq"]== 169)|(reference["name"]== "OG") & (reference["resname"]== "SER") & (reference["resseq"]== 169)|(reference["name"]== "OE1") & (reference["resname"]== "GLU") & (reference["resseq"]== 312)|(reference["name"]== "OE2") & (reference["resname"]== "GLU") & (reference["resseq"]== 312)|(reference["name"]== "NZ") & (reference["resname"]== "LYS") & (reference["resseq"]== 317)|(reference["name"]== "ND2") & (reference["resname"]== "ASN") & (reference["resseq"]== 350)|(reference["name"]== "OD1") & (reference["resname"]== "ASN") & (reference["resseq"]== 350)|(reference["name"]== "O") & (reference["resname"]== "ASN") & (reference["resseq"]== 350)|(reference["name"]== "OD1") & (reference["resname"]== "ASN") & (reference["resseq"]== 338)|(reference["name"]== "ND2") & (reference["resname"]== "ASN") & (reference["resseq"]== 338)|(reference["name"]== "O") & (reference["resname"]== "ASN") & (reference["resseq"]== 338)|(reference["name"]== "OD1") & (reference["resname"]== "ASN") & (reference["resseq"]== 155)|(reference["name"]== "ND2") & (reference["resname"]== "ASN") & (reference["resseq"]== 155)|(reference["name"]== "O") & (reference["resname"]== "ASN") & (reference["resseq"]== 155)|(reference["name"]== "O") & (reference["resname"]== "ASN") & (reference["resseq"]== 335)|(reference["name"]== "N") & (reference["resname"]== "ASN") & (reference["resseq"]== 335)|(reference["name"]== "HG1") & (reference["resname"]== "PRO") & (reference["resseq"]== 334)|(reference["name"]== "HG2") & (reference["resname"]== "PRO") & (reference["resseq"]== 334)|(reference["name"]== "HD1") & (reference["resname"]== "PRO") & (reference["resseq"]== 334)|(reference["name"]== "HD2") & (reference["resname"]== "PRO") & (reference["resseq"]== 334)|(reference["name"]== "O") & (reference["resname"]== "LEU") & (reference["resseq"]== 337)|(reference["name"]== "HZ1") & (reference["resname"]== "LYS") & (reference["resseq"]== 339)|(reference["name"]== "HZ2") & (reference["resname"]== "LYS") & (reference["resseq"]== 339)|(reference["name"]== "HZ3") & (reference["resname"]== "LYS") & (reference["resseq"]== 339)|(reference["name"]== "OD1") & (reference["resname"]== "ASP") & (reference["resseq"]== 96)|(reference["name"]== "OD2") & (reference["resname"]== "ASP") & (reference["resseq"]== 96)|(reference["name"]== "HA1") & (reference["resname"]== "GLY") & (reference["resseq"]== 333)|(reference["name"]== "HA2") & (reference["resname"]== "GLY") & (reference["resseq"]== 333)]
    filtered_ref.head()

    # cell 5

    filtered_ref.to_pickle(Path(file_stem) / 'xtal_waters_refposition.pkl')
    # cell 6

    colspecs_pdb = [(0, 6), (6, 11), (12, 16), (16, 17), (17, 20), (21, 22), (22, 26),
                (26, 27), (30, 38), (38, 46), (46, 54), (54, 60), (60, 66), (76, 78),
                (78, 80)]
    names_pdb = ['ATOM', 'serial', 'name', 'altloc', 'resname', 'chainid', 'resseq',
            'icode', 'x', 'y', 'z', 'occupancy', 'tempfactor', 'element', 'charge']

    pdb = pd.read_fwf(OW_pdb_path, names=names_pdb, colspecs=colspecs_pdb,skiprows=4)
    pdb
    #cell 7

    dist = 1.5  # Assume I can find a peak within 1.5 Angstrom of each chain G water
    for i, row in filtered_ref.iterrows():
        point1= [float(filtered_ref.loc[i,'x']),float(filtered_ref.loc[i,'y']),float(filtered_ref.loc[i,'z'])] 
        filtered_ref.loc[i,'OW_ID'] = 0 
        
        for q,raw in pdb.iterrows():
            point2= [float(pdb.loc[q,'x']),float(pdb.loc[q,'y']),float(pdb.loc[q,'z'])] 
            if (distance(point1,point2) < dist) :
                #replace position in reference dataframe 
                filtered_ref.loc[i,'x'] = pdb.loc[q,'x']
                filtered_ref.loc[i,'y'] = pdb.loc[q,'y']
                filtered_ref.loc[i,'z'] = pdb.loc[q,'z']
                filtered_ref.loc[i,'OW_ID'] = pdb.loc[q,'resseq'] 
                filtered_ref.loc[i,'OW_chain'] = pdb.loc[q,'chainid'] 

                break 

    #cell 8

    filtered_ref.to_pickle(Path(file_stem) / 'Waters_Odensity_xtalposition.pkl')
    finalgro = pd.read_pickle(Path(file_stem) / 'Waters_Odensity_xtalposition.pkl')

    # cell 9

    finalgro2 = finalgro.copy()
    x=0
    for i,row in finalgro2.iterrows():
            if x < 45:
                print('%s_%s_%s_%s' %(finalgro2.loc[i,'resname'],finalgro2.loc[i,'chainid'],finalgro2.loc[i,'name'], finalgro2.loc[i,'resseq']))
                finalgro2.loc[i,'resseq'] = '%s_%s_%s_%d' %(finalgro2.loc[i,'resname'],finalgro2.loc[i,'chainid'],finalgro2.loc[i,'name'], finalgro2.loc[i,'resseq'])
                x = x+1 
    #cell 10

    saveme = finalgro2[0:45][finalgro2[0:45]['resname']!='OEC']
    saveme.to_csv(Path(file_stem) / 'residues_used_for_hbond_analysis', sep='\t')

    #cell 11

    pdb_1 = pd.read_fwf(H1_pdb_path, names=names_pdb, colspecs=colspecs_pdb,skiprows=4)

    ##############
    pdb_2 = pd.read_fwf(H2_pdb_path, names=names_pdb, colspecs=colspecs_pdb,skiprows=4)

    # cell 12

    dist = 1.4
    for f, rof in finalgro2.iterrows():
        
        #average Ox pos
        point1= [float(finalgro2.loc[f,'x']),float(finalgro2.loc[f,'y']),float(finalgro2.loc[f,'z'])]  
        
        finalgro2.loc[f,'HW1_ID'] = 0
        finalgro2.loc[f,'HW1_x'] = 0
        finalgro2.loc[f,'HW1_y'] = 0
        finalgro2.loc[f,'HW1_z'] = 0
        
        finalgro2.loc[f,'HW2_ID'] = 0
        finalgro2.loc[f,'HW2_x'] = 0
        finalgro2.loc[f,'HW2_y'] = 0
        finalgro2.loc[f,'HW2_z'] = 0

    
        for g,rog in pdb_1.iterrows():
            point2= [pdb_1.loc[g,'x'],pdb_1.loc[g,'y'],pdb_1.loc[g,'z']] 
            if (distance(point1,point2) < dist) :
    #             print("HW1 Peak ID = %d " %pdb_1.loc[g,'resseq'])
                finalgro2.loc[f,'HW1_ID'] = pdb_1.loc[g,'resseq']
                finalgro2.loc[f,'HW1_x'] = pdb_1.loc[g,'x']
                finalgro2.loc[f,'HW1_y'] = pdb_1.loc[g,'y']
                finalgro2.loc[f,'HW1_z'] = pdb_1.loc[g,'z']
                break 
        for h,roh in pdb_2.iterrows():
            point3= [pdb_2.loc[h,'x'],pdb_2.loc[h,'y'],pdb_2.loc[h,'z']] 
            if (distance(point1,point3) < dist) :
                finalgro2.loc[f,'HW2_ID'] = pdb_2.loc[h,'resseq']
                finalgro2.loc[f,'HW2_x'] = pdb_2.loc[h,'x']
                finalgro2.loc[f,'HW2_y'] = pdb_2.loc[h,'y']
                finalgro2.loc[f,'HW2_z'] = pdb_2.loc[h,'z']

                break


    # cell 13

    finalgro2.to_pickle(Path(file_stem) / 'Waters_Odensity_xtalposition_withH.pkl')

    #cell 14

    finalgro2 = pd.read_pickle(Path(file_stem) / 'Waters_Odensity_xtalposition_withH.pkl')

    #Will make list of all matches with O-H distance < 3.0 Angstrom and can plot histogram for distribution 
    master_O_O_array = []
    master_O_H_array= []

    for p,rop in finalgro2.iterrows():
        
        ###################
        waterID_h1_partner_array =[]
        h1_partner_id_array =[]
        ow_h1_partner_dist_array=[]
        waterID_h1_partner_pdb_array=[]
        o_o_dist_h1_partner_array=[]
        
        waterID_h2_partner_array =[]
        h2_partner_id_array =[]
        ow_h2_partner_dist_array=[]
        waterID_h2_partner_pdb_array=[]
        o_o_dist_h2_partner_array=[]
        ###################
        
        dist_lowest = 3.0  #Max O-H distance 
        point1= [float(finalgro2.loc[p,'x']),float(finalgro2.loc[p,'y']),float(finalgro2.loc[p,'z'])]  #Ox

        for q,roq in finalgro2.iterrows(): 
            
            ###################Get H and Ox peaks for each element###################
            point2 =[float(finalgro2.loc[q,'HW1_x']),float(finalgro2.loc[q,'HW1_y']),float(finalgro2.loc[q,'HW1_z'])] 
            point2_ =[float(finalgro2.loc[q,'x']),float(finalgro2.loc[q,'y']),float(finalgro2.loc[q,'z'])] 
        
            if (distance(point1,point2) < dist_lowest) and (distance(point1, point2) < distance(point1, point2_)) and  (p!=q) and (finalgro2.loc[q,'OW_ID']!=0) and (finalgro2.loc[p,'OW_ID']!=0) :
                
                ###################If all conditions are met, append O-H1 bond instance to list###################
                waterID_h1_partner_array.append(finalgro2.loc[q,'resseq'])
                h1_partner_id_array.append(finalgro2.loc[q,'HW1_ID'])
                ow_h1_partner_dist_array.append(distance(point1,point2))                
                point4 =[float(finalgro2.loc[q,'x']),float(finalgro2.loc[q,'y']),float(finalgro2.loc[q,'z'])] #O
                o_o_dist_h1_partner_array.append(distance(point1,point4))
        
        ############ Insert array of all potential bonds which satisfy requirements ###################
        finalgro2.loc[p,'waterID_h1_partner']= str(waterID_h1_partner_array)
        finalgro2.loc[p,'h1_partner_id']= str(h1_partner_id_array)
        finalgro2.loc[p,'ow_h1_partner_dist'] =  str(ow_h1_partner_dist_array)
        finalgro2.loc[p,'o_o_dist_h1_partner'] = str(o_o_dist_h1_partner_array)

        #############################################################################################################
        
        for y,roy in finalgro2.iterrows(): 
            ###################Get H and Ox peaks for each element###################
            point3 =[float(finalgro2.loc[y,'HW2_x']),float(finalgro2.loc[y,'HW2_y']),float(finalgro2.loc[y,'HW2_z'])] 
            point3_ =[float(finalgro2.loc[y,'x']),float(finalgro2.loc[y,'y']),float(finalgro2.loc[y,'z'])] 
            
            ###################If all conditions are met, append O-H2 bond instance to list###################
            if (distance(point1,point3) < dist_lowest) and (distance(point1,point3) < distance(point1,point3_)) and (p!=y) and (finalgro2.loc[y,'OW_ID']!=0) and (finalgro2.loc[p,'OW_ID']!=0): 
                waterID_h2_partner_array.append(finalgro2.loc[y,'resseq'])
                h2_partner_id_array.append(finalgro2.loc[y,'HW2_ID'])
                ow_h2_partner_dist_array.append(distance(point1,point3))
                point5 =[float(finalgro2.loc[y,'x']),float(finalgro2.loc[y,'y']),float(finalgro2.loc[y,'z'])] 
                o_o_dist_h2_partner_array.append(distance(point1,point5))
        
        ############ Insert array of all potential bonds which satisfy requirements ###################        
        finalgro2.loc[p,'waterID_h2_partner']= str(waterID_h2_partner_array)
        finalgro2.loc[p,'h2_partner_id']= str(h2_partner_id_array)
        finalgro2.loc[p,'ow_h2_partner_dist'] =  str(ow_h2_partner_dist_array)
        finalgro2.loc[p,'o_o_dist_h2_partner'] = str(o_o_dist_h2_partner_array)
        
        ###################Append O-O and O-H distance info to master arrays for histogram plots ###################
        master_O_O_array.append(o_o_dist_h1_partner_array)
        master_O_O_array.append(o_o_dist_h2_partner_array)
        master_O_H_array.append(ow_h1_partner_dist_array)
        master_O_H_array.append(ow_h2_partner_dist_array)


    #cell 15

    O_O_flat = flatten_list(master_O_O_array)
    O_H_flat = flatten_list(master_O_H_array)

    plt.figure(figsize=[10,6])
    plt.hist(O_O_flat,bins=20, label='O-O distances')
    plt.hist(O_H_flat,bins=20, alpha=0.5, label='O-H distances')
    plt.legend(fontsize=20)
    plt.xlabel(r'Distance ($\AA$)')
    plt.xlim(1,4)
    plt.ylabel('Count')
    plt.savefig('Kanishk_MADI_distribution.png', dpi=300)
    print('***For true O-H distance***')
    print('33rd percentile = ', np.percentile(O_H_flat, 33))
    print('66th percentile = ', np.percentile(O_H_flat, 66))
    print('99th percentile = ', np.percentile(O_H_flat, 99))

    print('***For true O-O distance***')
    print('33rd percentile = ', np.percentile(O_O_flat, 33))
    print('66th percentile = ', np.percentile(O_O_flat, 66))
    print('99th percentile = ', np.percentile(O_O_flat, 99))
    # cell 16

    np.save(Path(file_stem) / 'O_O_master_array', O_O_flat)
    np.save(Path(file_stem) / 'O_H_master_array', O_H_flat)
    finalgro2.to_pickle(Path(file_stem) / 'hbond_analysis_3_20_22.pkl')

    # cell 17

    finalgro=pd.read_pickle(Path(file_stem) / 'hbond_analysis_3_20_22.pkl')
    finalgro=finalgro.drop(['HW1_x'], axis = 1)
    finalgro =finalgro.drop(['y'], axis = 1)
    finalgro = finalgro.drop(['z'], axis = 1)
    finalgro = finalgro.drop(['tempfactor'], axis = 1)
    finalgro = finalgro.drop(['occupancy'], axis = 1)
    finalgro = finalgro.drop(['HW1_y'], axis = 1)
    finalgro = finalgro.drop(['HW1_z'], axis = 1)
    finalgro = finalgro.drop(['HW2_x'], axis = 1)
    finalgro = finalgro.drop(['HW2_y'], axis = 1)
    finalgro = finalgro.drop(['HW2_z'], axis = 1)
    finalgro = finalgro.drop(['icode'], axis = 1)
    finalgro = finalgro.drop(['altloc'], axis = 1)
    finalgro = finalgro.drop(['element'], axis = 1)
    finalgro = finalgro.drop(['charge'], axis = 1)
    finalgro = finalgro.replace('NaN', 0)
    finalgro = finalgro.replace(np.nan, 0)

    # cell 18

    ############Initialize columns which will contain final O-H information############
    finalgro['ow_h1_partner_dist_fin'] = 0
    finalgro['ow_h2_partner_dist_fin'] = 0
    finalgro['waterID_h1_partner_fin'] = 0
    finalgro['waterID_h2_partner_fin'] = 0
    finalgro['o_o_dist_h1_partner_fin'] = 0
    finalgro['o_o_dist_h2_partner_fin'] = 0



    for f, rof in finalgro.iterrows():
        
        ############Grab all info about first element############
        res1 = finalgro.loc[f,'resseq']

        case1 = finalgro.loc[f,'chainid'].isupper()
        h1_dist1 = finalgro.loc[f,'ow_h1_partner_dist']
        h1_dist1 = h1_dist1.replace('[', '')
        h1_dist1 = h1_dist1.replace(']', '')
        h1_dist1  = h1_dist1.replace(" ", "")
        h1_dist1 = h1_dist1.split(',')
        
        h1_1 = finalgro.loc[f,'waterID_h1_partner']
        h1_1 = h1_1.replace('[', '')
        h1_1 = h1_1.replace(']', '')
        h1_1  = h1_1.replace(" ", "")
        h1_1 = h1_1.split(',')
        
        h2_1=finalgro.loc[f,'waterID_h2_partner']
        h2_1 = h2_1.replace('[', '')
        h2_1 = h2_1.replace(']', '')
        h2_1  = h2_1.replace(" ", "")
        h2_1 = h2_1.split(',')

        h2_dist1 = finalgro.loc[f,'ow_h2_partner_dist']
        h2_dist1 = h2_dist1.replace('[', '')
        h2_dist1 = h2_dist1.replace(']', '')
        h2_dist1  = h2_dist1.replace(" ", "")
        h2_dist1 = h2_dist1.split(',')
        
        O_Odist1 = finalgro.loc[f,'o_o_dist_h1_partner']
        O_Odist1 = O_Odist1.replace('[', '')
        O_Odist1 = O_Odist1.replace(']', '')
        O_Odist1  = O_Odist1.replace(" ", "")
        O_Odist1 = O_Odist1.split(',')
        
        O_Odist2 = finalgro.loc[f,'o_o_dist_h2_partner']
        O_Odist2 = O_Odist2.replace('[', '')
        O_Odist2 = O_Odist2.replace(']', '')
        O_Odist2 = O_Odist2.replace(" ", "")
        O_Odist2 = O_Odist2.split(',')
        

        ############################################################

        #############Convert O-H distances to array of floats and sort ##################
        try:
            h1_dist1 = [float(i) for i in h1_dist1]
        except ValueError:
            h1_dist1 =[]
        try:
            h2_dist1 = [float(i) for i in h2_dist1]
        except ValueError:
            h2_dist1 =[]
            
        h1_1 = [x for y, x in sorted(zip(h1_dist1, h1_1))]
        h2_1 = [x for _,x in sorted(zip(h2_dist1,h2_1))]
            
        O_Odist1 = [float(x) for _,x in sorted(zip(h1_dist1,O_Odist1))]
        O_Odist2 = [float(x) for _,x in sorted(zip(h2_dist1,O_Odist2))]
        
        h1_dist1 = [float(x) for _,x in sorted(zip(h1_dist1,h1_dist1))]
        h2_dist1 = [float(x) for _,x in sorted(zip(h2_dist1,h2_dist1))]
            
    ############Grab all info about second element############
        for z, rox in finalgro.iterrows():
        
            case2 =   finalgro.loc[z,'chainid'].isupper()
            
            if case1 !=case2: #If comparing two separate chains, skip. Only compare within same chain. 
                continue
                
            res2 = finalgro.loc[z,'resseq']

            h1_dist2 = finalgro.loc[z,'ow_h1_partner_dist']
            h1_dist2 = h1_dist2.replace('[', "")
            h1_dist2 = h1_dist2.replace(']', "")
            h1_dist2 = h1_dist2.replace(" ", "")
            h1_dist2 = h1_dist2.split(',')
        
            h2_dist2 = finalgro.loc[z,'ow_h2_partner_dist']
            h2_dist2 = h2_dist2.replace('[', "")
            h2_dist2 = h2_dist2.replace(']', "")
            h2_dist2 = h2_dist2.replace(" ", "")
            h2_dist2 = h2_dist2.split(',')
            
            h1_2 = finalgro.loc[z,'waterID_h1_partner']
            h1_2 = h1_2.replace('[', "")
            h1_2 = h1_2.replace(']', "")
            h1_2 = h1_2.replace(" ", "")
            h1_2 = h1_2.split(',')
            
            h2_2=finalgro.loc[z,'waterID_h2_partner']
            h2_2 = h2_2.replace('[', "")
            h2_2 = h2_2.replace(']', "")
            h2_2 = h2_2.replace(" ", "")
            h2_2 = h2_2.split(',')
        

            #############Convert O-H distances to array of floats and sort ##################
            try:
                h1_dist2 = [float(i) for i in h1_dist2]
            except ValueError:  #This error generally means that you have an empty array 
                h1_dist2=[]
            try:
                h2_dist2 = [float(i) for i in h2_dist2]
            except ValueError:
                h2_dist2=[] 
                    
            h1_2 = [x for y, x in sorted(zip(h1_dist2, h1_2))]
            h2_2 = [x for _,x in sorted(zip(h2_dist2,h2_2))]
        
            h1_dist2 = [float(x) for _,x in sorted(zip(h1_dist2,h1_dist2))]
            h2_dist2 = [float(x) for _,x in sorted(zip(h2_dist2,h2_dist2))]
            
            
            #############Take strongest bond from each case and see if there are duplicates#############
            ############# If there are duplicates, then find *true* direction of O-H bond  #############
            try:
                h1_dist1_ = h1_dist1[0]
                h1_1_=h1_1[0]
            except IndexError:
                h1_1_=0
                h1_dist1_  = 1000000
            
            try:
                h2_dist1_ = h2_dist1[0]
                h2_1_ = h2_1[0]
            except IndexError:
                h2_dist1_ = 1000000
                h2_1_ = 0
            
            try:
                h1_dist2_ = h1_dist2[0]
                h1_2_=h1_2[0]
            except IndexError:
                h1_dist2_ = 1000000  
                h1_2_ = 0
            
            try:
                h2_dist2_ = h2_dist2[0]
                h2_2_=h2_2[0]
            except IndexError:
                h2_dist2_ = 1000000 
                h2_2_ = 0
            
            try:
                h1_1_ = h1_1_.replace("'", "")
            except AttributeError:
                pass
            try:
                h2_1_ = h2_1_.replace("'", "")
            except AttributeError:
                pass
            try:
                h1_2_ = h1_2_.replace("'", "")
            except AttributeError:
                pass
                
            try:
                h2_2_ = h2_2_.replace("'", "")
            except AttributeError:
                pass
                
            if (str(res1) == str(h1_2_)) and (str(res2) == str(h1_1_)):

                distarray =[h1_dist1_, h1_dist2_]
                if h1_dist1_ == np.min(distarray): 
                    try:
                        h1_2 = h1_2[1:]
                        h1_dist2 = h1_dist2[1:]

                    except IndexError : 
                        h1_2 =[0]
                        h1_dist2 = [0]   


                elif h1_dist2_ == np.min(distarray): 
                    try:
                        h1_1 = h1_1[1:]
                        h1_dist1 = h1_dist1[1:]                    
                    
                    except IndexError : 
                        h1_1 = [0]
                        h1_dist1 = [0]                  

    
            if (str(res1) == str(h2_2_)) and (str(res2) == str(h2_1_)):
                distarray =[h2_dist1_, h2_dist2_]
                if h2_dist1_ == np.min(distarray): 
                    try:
                        h2_2 = h2_2[1:]
                        h2_dist2 = h2_dist2[1:]    
                        
                    except IndexError:
                        h2_2 = [0]
                        h2_dist2 = [0]   

                elif h2_dist2_ == np.min(distarray): 
                    try:
                        h2_1 = h2_1[1:]
                        h2_dist1 = h2_dist1[1:]                   
                    
                    except IndexError: 
                        h2_1 = [0]
                        h2_dist1 = [0]                    

            if (str(res1) == str(h2_2_)) and (str(res2) == str(h1_1_)):

                distarray =[h1_dist1_, h2_dist2_]
                if h1_dist1_ == np.min(distarray): 
                    try:
                        h2_2 = h2_2[1:]
                        h2_dist2 = h2_dist2[1:] 

                    except IndexError: 
                        h2_2 = [0]
                        h2_dist2= [0]       

                    print('chain' + str(finalgro.loc[z,'chainid']) + ' ' + str(res1) + ' and  ' + str (h1_1_) + ' (h1) beat ' + str(res2) + ' and ' + str (h2_2_)+ '(h2)')
            
                elif h2_dist2_ == np.min(distarray): 
                    try:
                        h1_1 = h1_1[1:]
                        h1_dist1 = h1_dist1[1:]                                   
                    except IndexError: 
                        h1_1 = [0]
                        h1_dist1 = [0]                   
                    print('chain' + str(finalgro.loc[z,'chainid']) + ' ' + str(res2) + ' and  ' + str (h2_2_) + ' (h2) beat ' + str(res1) + ' and ' + str (h1_1_)+ '(h1)')
                        
                
            if (str(res1) == str(h1_2_)) and (str(res2) == str(h2_1_)):

                distarray =[h2_dist1_, h1_dist2_]
                if h2_dist1_ == np.min(distarray): 
                    try:
                        h1_dist2 = h1_2[1:]
                        h1_dist2 = h1_dist2[1:]  
                    except IndexError: 
                        h1_dist2 = [0]
                        h1_dist2 =  [0] 

                        
                    print('chain' + str(finalgro.loc[z,'chainid']) +  ' ' + str(res1) + ' and  ' + str (h2_1_) + ' (h2) beat '+ str(res2) + ' and ' + str (h1_2_)+ '(h1)')
            
                elif h1_dist2_ == np.min(distarray): 
                    try:
                        h2_1 = h2_1[1:]
                        h2_dist1 = h2_dist1[1:]                                  
                    except IndexError: 
                        h2_1 = [0]
                        h2_dist1 = [0]                    

                    print('chain' + str(finalgro.loc[z,'chainid']) + ' ' + str(res2) + ' and  ' + str (h1_2_) + ' (h1) beat ' +str(res1) + ' and ' + str (h2_1_)+ '(h2)')
        
            ########################Finalizing O-H bonds########################
        
        if ((len(h2_1) ==0) and (len(h1_1)>0)):
        
            finalgro.loc[f,'waterID_h2_partner_fin'] = 0
            finalgro.loc[f,'ow_h2_partner_dist_fin'] = 0
            finalgro.loc[f,'o_o_dist_h2_partner_fin'] = 0

            finalgro.loc[f,'waterID_h1_partner_fin'] = [h1_1[0]]
            finalgro.loc[f,'ow_h1_partner_dist_fin'] = [h1_dist1[0]]
            finalgro.loc[f,'o_o_dist_h1_partner_fin'] = O_Odist1[0]


        elif ((len(h1_1) ==0) and (len(h2_1)>0)):

            finalgro.loc[f,'waterID_h2_partner_fin'] = [h2_1[0]]
            finalgro.loc[f,'ow_h2_partner_dist_fin'] = [h2_dist1[0]]
            finalgro.loc[f,'o_o_dist_h2_partner_fin'] = O_Odist2[0]

            finalgro.loc[f,'waterID_h1_partner_fin'] = 0
            finalgro.loc[f,'ow_h1_partner_dist_fin'] = 0
            finalgro.loc[f,'o_o_dist_h1_partner_fin'] =0
        ######If same ID is not registered for H1 and H2 , this is easy ######

        elif ( (len(h1_1) >0) and (len(h2_1)>0) and ((h2_1[0] != h1_1[0]) )):

            finalgro.loc[f,'waterID_h2_partner_fin'] = [h2_1[0]]
            finalgro.loc[f,'ow_h2_partner_dist_fin'] = [h2_dist1[0]]
            finalgro.loc[f,'o_o_dist_h2_partner_fin'] = O_Odist2[0]

            finalgro.loc[f,'waterID_h1_partner_fin'] = [h1_1[0]]
            finalgro.loc[f,'ow_h1_partner_dist_fin'] = [h1_dist1[0]]
            finalgro.loc[f,'o_o_dist_h1_partner_fin'] = O_Odist1[0]
            
        ######However, if H1 and H2 are same, do some checks######  

        else :  

            #If you have backup o-h2 bonds 
            if ((len(h1_dist1) == 1) and  (len(h1_dist1) < len(h2_dist1))):

                finalgro.loc[f,'waterID_h2_partner_fin'] = [h2_1[1]]
                finalgro.loc[f,'ow_h2_partner_dist_fin'] = [h2_dist1[1]]
                finalgro.loc[f,'o_o_dist_h2_partner_fin'] = O_Odist2[1]

                finalgro.loc[f,'waterID_h1_partner_fin'] = [h1_1[0]]
                finalgro.loc[f,'ow_h1_partner_dist_fin'] = [h1_dist1[0]]
                finalgro.loc[f,'o_o_dist_h1_partner_fin'] = O_Odist1[0]

            #If you have backup o-h1 bonds          
            elif ((len(h2_dist1) == 1) and len(h2_dist1) < len(h1_dist1)): 

                finalgro.loc[f,'waterID_h2_partner_fin'] = [h2_1[0]]
                finalgro.loc[f,'ow_h2_partner_dist_fin'] = [h2_dist1[0]]
                finalgro.loc[f,'o_o_dist_h2_partner_fin'] = O_Odist2[0]

                finalgro.loc[f,'waterID_h1_partner_fin'] = [h1_1[1]]
                finalgro.loc[f,'ow_h1_partner_dist_fin'] = [h1_dist1[1]]
                finalgro.loc[f,'o_o_dist_h1_partner_fin'] = O_Odist1[1]

            elif ( (len(h1_1) >0) and (len(h2_1)>0) and (h2_1[0] == h1_1[0])): 

                if h1_dist1[0] < h2_dist1[0]: 
                    finalgro.loc[f,'waterID_h1_partner_fin'] = [h1_1[0]]
                    finalgro.loc[f,'ow_h1_partner_dist_fin'] = [h1_dist1[0]]
                    finalgro.loc[f,'o_o_dist_h1_partner_fin'] = O_Odist1[0]
                    try:
                        finalgro.loc[f,'waterID_h2_partner_fin'] = [h2_1[1]]
                        finalgro.loc[f,'ow_h2_partner_dist_fin'] = [h2_dist1[1]]
                        finalgro.loc[f,'o_o_dist_h2_partner_fin'] = O_Odist2[1]
                            
                    except IndexError:

                        finalgro.loc[f,'waterID_h2_partner_fin'] = 0
                        finalgro.loc[f,'ow_h2_partner_dist_fin'] = 0
                        finalgro.loc[f,'o_o_dist_h2_partner_fin'] = 0

                elif h2_dist1[0] < h1_dist1[0]: 
                    finalgro.loc[f,'waterID_h2_partner_fin'] = [h2_1[0]]
                    finalgro.loc[f,'ow_h2_partner_dist_fin'] = [h2_dist1[0]]
                    finalgro.loc[f,'o_o_dist_h2_partner_fin'] = O_Odist2[0]
                    try:
                        finalgro.loc[f,'waterID_h1_partner_fin'] = [h1_1[1]]
                        finalgro.loc[f,'ow_h1_partner_dist_fin'] = [h1_dist1[1]]
                        finalgro.loc[f,'o_o_dist_h1_partner_fin'] = O_Odist1[0]

                    except IndexError: #If there are no backup options 
                        finalgro.loc[f,'waterID_h1_partner_fin'] = 0
                        finalgro.loc[f,'ow_h1_partner_dist_fin'] =0
                        finalgro.loc[f,'o_o_dist_h1_partner_fin'] = 0

    # cell 19

    for f, rof in finalgro.iterrows():
            
        ############Grab all info about first element############
        res1 = finalgro.loc[f,'resseq']
        case1 = finalgro.loc[f,'chainid'].isupper()
        h1_dist1 = finalgro.loc[f,'ow_h1_partner_dist_fin']
        h2_dist1 = finalgro.loc[f,'ow_h2_partner_dist_fin']
        
        h1_1 = finalgro.loc[f,'waterID_h1_partner_fin']
        h1_1 = str(h1_1)

        h2_1=finalgro.loc[f,'waterID_h2_partner_fin']
        h2_1=str(h2_1)

        O_Odist1 = finalgro.loc[f,'o_o_dist_h1_partner_fin']


        O_Odist2 = finalgro.loc[f,'o_o_dist_h2_partner_fin']

    #  ############Grab all info about second element############
        for z, rox in finalgro.iterrows():

            case2 =   finalgro.loc[z,'chainid'].isupper()

            if case1 !=case2: #If comparing two separate chains, skip. Only compare within same chain. 
                continue

            res2 = finalgro.loc[z,'resseq']

            h1_dist2 = finalgro.loc[z,'ow_h1_partner_dist_fin']
            h2_dist2 = finalgro.loc[z,'ow_h2_partner_dist_fin']
            h1_2 = finalgro.loc[z,'waterID_h1_partner_fin']
            h2_2=finalgro.loc[z,'waterID_h2_partner_fin']
            h1_2 = str(h1_2)
            h2_2 = str(h2_2)
            
            h1_1 = h1_1.replace("'", "")
            h2_1 = h2_1.replace("'", "")
            h1_2 = h1_2.replace("'", "")
            h2_2 = h2_2.replace("'", "")


            ############# If there are duplicates, then find *true* direction of O-H bond  #############

            if (str(res1) == str(h1_2)) and (str(res2) == str(h1_1)):
        
                distarray =[h1_dist1, h1_dist2]
                if h1_dist1 == np.min(distarray): 

                    finalgro.loc[z,'waterID_h1_partner_fin'] =0          

                elif h1_dist2 == np.min(distarray): 

                    finalgro.loc[f,'waterID_h1_partner_fin'] =0  

            if (str(res1) == str(h2_2)) and (str(res2) == str(h2_1)):
    
                distarray =[h2_dist1, h2_dist2]
                if h2_dist1 == np.min(distarray): 

                    finalgro.loc[z,'waterID_h2_partner_fin'] =0            

                elif h2_dist2 == np.min(distarray): 

                    finalgro.loc[f,'waterID_h2_partner_fin'] =0  
            if (str(res1) == str(h2_2)) and (str(res2) == str(h1_1)):
    
                distarray =[h1_dist1, h2_dist2]
                if h1_dist1 == np.min(distarray): 

                    finalgro.loc[z,'waterID_h2_partner_fin'] =0  

                elif h2_dist2 == np.min(distarray): 

                    finalgro.loc[f,'waterID_h1_partner_fin'] =0 



            if (str(res1) == str(h1_2)) and (str(res2) == str(h2_1)):

                distarray =[h2_dist1, h1_dist2]
                if h2_dist1 == np.min(distarray): 

                    finalgro.loc[z,'waterID_h1_partner_fin'] =0 


                elif h1_dist2 == np.min(distarray): 
                    finalgro.loc[f,'waterID_h2_partner_fin'] =0 
    # cell 20

    finalgro.to_pickle(Path(file_stem) / 'hbond_analysis_filtered.pkl')

    # cell 21

    finalgro=pd.read_pickle(Path(file_stem) / 'hbond_analysis_filtered.pkl')
    finalgro = finalgro.replace('NaN', 0)
    finalgro = finalgro.replace(np.nan, 0)
    finalgro.head()

    # cell 22

    from collections import Counter

    a = np.array(finalgro[0:152]['waterID_h1_partner_fin'].values)
    a_oh = np.array(finalgro[0:152]['ow_h1_partner_dist_fin'].values)
    b = np.array(finalgro[0:152]['waterID_h2_partner_fin'].values)
    b_oh = np.array(finalgro[0:152]['ow_h2_partner_dist_fin'].values)

    master_array =[]
    o_h_distances=[]
    h_type=[]


    for x in range(0,len(a)):
        master_array.append(a[x])
        o_h_distances.append(a_oh[x])
        h_type.append(1)
        
    for y in range(0,len(b)):
        master_array.append(b[y])
        o_h_distances.append(b_oh[y])
        h_type.append(2)


    for x in range(0,len(master_array)):
        dist_array=[]
        h_array=[]
        index_array=[]
        if master_array[x] ==0:
            continue 
        count = 1 
        dist_array.append(o_h_distances[x])
        h_array.append(h_type[x])
        index_array.append(x)
        for i in range(0,len(master_array)):
        
            try:
                if (float(master_array[x]) == float(master_array[i])) and (i!=x): 
                    count = count +1
                    dist_array.append(o_h_distances[i])
                    h_array.append(h_type[i])
                    index_array.append(i)
            except (TypeError, ValueError): 
                if ((master_array[x]) == (master_array[i])) and (i!=x) : 
                    count = count +1
                    dist_array.append(o_h_distances[i])
                    h_array.append(h_type[i])
                    index_array.append(i)

        if count > 2:
            print('%d counts of %s' %(count, master_array[x]))

            dropindex_location= np.array(dist_array).argmax()
            dropindex=index_array[dropindex_location]
            drophtype = h_array[dropindex_location]

            if drophtype ==1:
                finalgro['waterID_h1_partner_fin'].iloc[dropindex]=0
            else:
                finalgro['waterID_h2_partner_fin'].iloc[dropindex-152]=0
                

    # cell 23

    #This function is used when you have a double acceptor situation 
    def plot_dub(z, p,chain):
        strength.append(finalgro['ow_h1_partner_dist_fin'].iloc[p])
        ax1.scatter(z,finalgro['ow_h1_partner_dist_fin'].iloc[p],s=100, color='blue',edgecolor='black')
        ax1.scatter(z,finalgro['o_o_dist_h1_partner_fin'].iloc[p] ,s=100,color='orange',edgecolor='black')
        ax2.plot(z, 1,color='white') 
        ticklabels1.append(finalgro['resseq'].iloc[p])
        ticklabels2.append(finalgro['waterID_h1_partner_fin'].iloc[p])
        prev1 = finalgro['waterID_h1_partner_fin'].iloc[p]
        z=z+1
        strength.append(finalgro['ow_h2_partner_dist_fin'].iloc[p])
        ax1.scatter(z,finalgro['ow_h2_partner_dist_fin'].iloc[p], s=100, edgecolor='black', color='blue')
        print('Double acceptor:', finalgro['resseq'].iloc[p])
        print('Partner IDs:',finalgro['waterID_h1_partner_fin'].iloc[p] ,  finalgro['waterID_h2_partner_fin'].iloc[p])
        print('***')
        ax1.scatter(z,finalgro['o_o_dist_h2_partner_fin'].iloc[p] ,s=100, edgecolor='black',color='orange')
        ax2.plot(z, 1,color='white') 
        ticklabels1.append(finalgro['resseq'].iloc[p])
        ticklabels2.append(finalgro['waterID_h2_partner_fin'].iloc[p])
        prev1 = finalgro['waterID_h2_partner_fin'].iloc[p]
        return z,prev1

    # cell 24

    chain = 'G'
    chan = 'Cl'

    # cell 25

    donor=[]
    acceptor =[]
    strength =[]

    OW_peaks_involved = []
    H1_peaks_involved = []
    H2_peaks_involved = []

    H1_partners_involved = []
    H2_partners_involved = []

    fig = plt.figure(figsize=(40, 5))
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Water ID (acceptor)')
    ax2 = ax1.twiny()
    ax2.set_xlabel('Water ID (donor)' )
    ax2.xaxis.labelpad=10.0
    ax1.set_ylabel(r'Distance $(\AA$)')

    ticklabels2=[]
    ticklabels1=[]
    strength = []

    st_O1 =  ['ASP_A_O_61', 'ASP_A_CB_61', 'ASP_A_OD1_61', 'ASP_A_OD2_61', 'GLU_A_OE1_65' , 'GLU_A_OE2_65', 'SER_A_O_169', 'SER_A_OG_169', 'GLY_A_O_171','ASN_A_O_181','ASN_A_OD1_181','ASN_A_ND2_181','ARG_A_NH1_334','ARG_A_NH2_334','ASN_A_N_335','ASN_A_O_335','ASN_A_O_338','ASN_A_OD1_338','ASN_A_ND2_338','PRO_A_O_340','OEC_A_O1_601','OEC_A_O2_601','OEC_A_O3_601','OEC_A_O4_601','OEC_A_O5_601'           ,'OEC_A_W1_601','OEC_A_W2_601','OEC_A_W3_601','OEC_A_W4_601','SER_B_O_169','SER_B_OG_169','GLY_B_CA_333','ASN_C_O_155','ASN_C_OD1_155','ASN_C_ND2_155','GLY_C_O_171','GLY_C_CA_333','PRO_C_CG_334','PRO_C_CD_334','THR_C_OG1_335','THR_C_CG2_335','LEU_C_O_337','LYS_C_NZ_339','GLU_D_CG_312','GLU_D_OE1_312','GLU_D_OE2_312','LYS_D_NZ_317','ASN_D_O_338','ASN_D_OD1_338','ASN_D_ND2_338','ASN_D_O_350','ASN_D_OD1_350','ASN_D_ND2_350', 26,27,28, 29, 30, 31, 32, 33,  34, 35,36,37,38, 39, 76,77, 101, 102, 103, 104, 105, 106, 107, 108]

    st_Cl = [  'ASP_A_O_61', 'ASP_A_CB_61', 'ASP_A_OD1_61', 'ASP_A_OD2_61', 'GLU_A_OE1_65' , 'GLU_A_OE2_65', 'SER_A_O_169', 'SER_A_OG_169', 'GLY_A_O_171','ASN_A_O_181','ASN_A_OD1_181','ASN_A_ND2_181','ARG_A_NH1_334','ARG_A_NH2_334','ASN_A_N_335','ASN_A_O_335','ASN_A_O_338','ASN_A_OD1_338','ASN_A_ND2_338','PRO_A_O_340','OEC_A_O1_601','OEC_A_O2_601','OEC_A_O3_601','OEC_A_O4_601','OEC_A_O5_601'           ,'OEC_A_W1_601','OEC_A_W2_601','OEC_A_W3_601','OEC_A_W4_601','SER_B_O_169','SER_B_OG_169','GLY_B_CA_333','ASN_C_O_155','ASN_C_OD1_155','ASN_C_ND2_155','GLY_C_O_171','GLY_C_CA_333','PRO_C_CG_334','PRO_C_CD_334','THR_C_OG1_335','THR_C_CG2_335','LEU_C_O_337','LYS_C_NZ_339','GLU_D_CG_312','GLU_D_OE1_312','GLU_D_OE2_312','LYS_D_NZ_317','ASN_D_O_338','ASN_D_OD1_338','ASN_D_ND2_338','ASN_D_O_350','ASN_D_OD1_350','ASN_D_ND2_350', 59,60,66,67,68,69,70,21, 22,23,24,25,61,62,  40,  41, 42,  150,   117, 119, 121,122,150]

    st_O4 =[ 'ASP_A_O_61', 'ASP_A_CB_61', 'ASP_A_OD1_61', 'ASP_A_OD2_61', 'SER_A_O_169', 'SER_A_OG_169', 'GLY_A_O_171','ASN_A_O_181','ASN_A_OD1_181','ASN_A_ND2_181','ARG_A_NH1_334','ARG_A_NH2_334','ASN_A_N_335','ASN_A_O_335','ASN_A_O_338','ASN_A_OD1_338','ASN_A_ND2_338','PRO_A_O_340','OEC_A_O1_601','OEC_A_O2_601','OEC_A_O3_601','OEC_A_O4_601','OEC_A_O5_601'           ,'OEC_A_W1_601','OEC_A_W2_601','OEC_A_W3_601','OEC_A_W4_601','SER_B_O_169','SER_B_OG_169','GLY_B_CA_333','ASN_C_O_155','ASN_C_OD1_155','ASN_C_ND2_155','GLY_C_O_171','GLY_C_CA_333','PRO_C_CG_334','PRO_C_CD_334','THR_C_OG1_335','THR_C_CG2_335','LEU_C_O_337','LYS_C_NZ_339','GLU_D_CG_312','GLU_D_OE1_312','GLU_D_OE2_312','LYS_D_NZ_317','ASN_D_O_338','ASN_D_OD1_338','ASN_D_ND2_338','ASN_D_O_350','ASN_D_OD1_350','ASN_D_ND2_350', 19,20,48,49, 50, 51, 52, 53, 71, 72, 73]
    
    print("done")
    donor, acceptor, strength = [], [], []
    ticklabels1, ticklabels2 = [], []
    OW_peaks_involved, H1_peaks_involved, H2_peaks_involved = [], [], []
    H1_partners_involved, H2_partners_involved = [], []
    
    O1(chain, finalgro, st_O1, ax1, ax2, donor, acceptor, strength, ticklabels1, ticklabels2, OW_peaks_involved, H1_peaks_involved, H2_peaks_involved, H1_partners_involved, H2_partners_involved, plot_dub, chan = 'O1')
    
    # 2. Reset with fresh lists for Cl
    donor, acceptor, strength = [], [], []
    ticklabels1, ticklabels2 = [], []
    OW_peaks_involved, H1_peaks_involved, H2_peaks_involved = [], [], []
    H1_partners_involved, H2_partners_involved = [], []

    Cl(chain, finalgro, st_Cl, ax1, ax2, donor, acceptor, strength, ticklabels1, ticklabels2, OW_peaks_involved, H1_peaks_involved, H2_peaks_involved, H1_partners_involved, H2_partners_involved, plot_dub, chan = 'Cl')
    
    # 3. Reset with fresh lists for O4
    donor, acceptor, strength = [], [], []
    ticklabels1, ticklabels2 = [], []
    OW_peaks_involved, H1_peaks_involved, H2_peaks_involved = [], [], []
    H1_partners_involved, H2_partners_involved = [], []

    O4(chain, finalgro, st_O4, ax1, ax2, donor, acceptor, strength, ticklabels1, ticklabels2, OW_peaks_involved, H1_peaks_involved, H2_peaks_involved, H1_partners_involved, H2_partners_involved, plot_dub, chan = 'O4')
    
    print("even more done!")

def pymol_O4():
    ticklabels2 = np.load(Path(pickle_stem) / 'ticklabels2_O4.npy')
    ticklabels1 = np.load(Path(pickle_stem) / 'ticklabels1_O4.npy')
    ticklabels = [i.strip("'") for i in ticklabels1] + [i.strip("'") for i in ticklabels2]

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
    waters += [50,51,52,53,73]

    
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
        yellow_waters = np.load(Path(pickle_stem) / "yellow_water_indices.npy")
        existing_selections = cmd.get_names("selections")

        #------COLOR YELLOW--------
        for water in yellow_waters: 
            sel_name = f"water_{water}"
            
            if sel_name in existing_selections:
                cmd.color("yellow", sel_name)
                cmd.set("label_color", "black", sel_name)
                

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
        arrows = np.load(Path(pickle_stem) / 'arrows_O4.npy')

        for i, (start, end, color) in enumerate(arrows):
            cgo_arrow(start, end, color=color, radius=0.2, hlength=0.4, hradius=0.4, gap=0.6, name=f'arrow_{i+1}')
        #----PEAK PDB FILES-----
        cmd.load(H1_peaks)
        cmd.load(H2_peaks)
        cmd.load(OW_peaks)

        #----SHOW SPHERES-----
        spheres = np.load(Path(pickle_stem) / 'spheres_O4.npy')
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

def pymol_Cl():
    ticklabels2 = np.load(Path(pickle_stem) / 'ticklabels2_Cl.npy')
    ticklabels1 = np.load(Path(pickle_stem) / 'ticklabels1_Cl.npy')
    ticklabels = [i.strip("'") for i in ticklabels1] + [i.strip("'") for i in ticklabels2]

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
    waters += [25,67,68,69,70,126,121,122,117,119,129]
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
                cmd.set("label_color", "black", sel_name)
                

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
        arrows = np.load(Path(pickle_stem) / "arrows_Cl.npy")

        for i, (start, end, color) in enumerate(arrows):
            cgo_arrow(start, end, color=color, radius=0.2, hlength=0.4, hradius=0.4, gap=0.6, name=f'arrow_{i+1}')
        #----PEAK PDB FILES-----
        cmd.load(H1_peaks)
        cmd.load(H2_peaks)
        cmd.load(OW_peaks)

        #----SHOW SPHERES-----
        spheres = np.load(Path(pickle_stem) / 'spheres_Cl.npy')
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

def pymol_O1():
    ticklabels2 = np.load(Path(pickle_stem) / 'ticklabels2_O1.npy')
    ticklabels1 = np.load(Path(pickle_stem) / 'ticklabels1_O1.npy')
    print(ticklabels1)
    print(ticklabels2)
    ticklabels = [i.strip("'") for i in ticklabels1] + [i.strip("'") for i in ticklabels2]

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
    waters += [106,107,35]
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
        yellow_waters = np.load(Path(pickle_stem) / "yellow_water_indices.npy")
        existing_selections = cmd.get_names("selections")

        #------COLOR YELLOW--------
        for water in yellow_waters: 
            sel_name = f"water_{water}"
            
            if sel_name in existing_selections:
                cmd.color("yellow", sel_name)
                cmd.set("label_color", "black", sel_name)
                

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
        arrows = np.load(Path(pickle_stem) / "arrows_O1.npy")

        for i, (start, end, color) in enumerate(arrows):
            cgo_arrow(start, end, color=color, radius=0.2, hlength=0.4, hradius=0.4, gap=0.6, name=f'arrow_{i+1}')
        #----PEAK PDB FILES-----
        cmd.load(H1_peaks)
        cmd.load(H2_peaks)
        cmd.load(OW_peaks)

        #----SHOW SPHERES-----
        spheres = np.load(Path(pickle_stem) / "spheres_O1.npy")

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

def main():
    try:
        import pymol
        from pymol import cmd, cgo, CmdException
        IN_PYMOL = True
    except ImportError:
        IN_PYMOL = False

    if IN_PYMOL:
        print("Running inside PyMOL")
        pymol_Cl()


    else:
        print("Running from a standard terminal")
        MADI_analysis()


main()


""" Detects if run in pymol vs terminal, maybe consider using?
import sys

try:
    import pymol
    from pymol import cmd, cgo, CmdException
    IN_PYMOL = True
except ImportError:
    IN_PYMOL = False

if IN_PYMOL:
    print("Running inside PyMOL")

else:
    print("Running from a standard terminal")

"""