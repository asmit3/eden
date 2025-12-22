import numpy as np
import matplotlib.pyplot as plt
import sys

fig, ax = plt.subplots()

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)
pair=sys.argv[1]
stepsize=0.01
timepoints = ['0F','0F_L98', '2F','2F_L98']
steps = np.arange(-1, 2, stepsize)
all_map_values_early = []
all_map_values_mid = []
linewidth = 3
for timepoint in timepoints:
  map_values = []
  dist_values = []
  with open('%s_lm_%s_map_values_noNormalization.csv'%(timepoint,pair), 'r') as fin:
    for line in fin:
      if line !='\n':
        bx = line.split(',')
        if 'atoms' in line:
          atom1=bx[1].strip()
          atom2=bx[2].strip()
          continue
        map_values.append(float(bx[1]))
        dist_values.append(float(bx[0]))

    map_values.reverse()
    dist_values.reverse()
    dist_values = [-x for x in dist_values]
    #if timepoint =='2F':
    if timepoint == '0F':
      ax.plot(dist_values, [x+0 for x in map_values],'--', label=timepoint, c='blue', linewidth=linewidth)
    elif timepoint == '0F_L98':
      ax.plot(dist_values, [x+0 for x in map_values], label=timepoint, c='blue', linewidth=linewidth)
    elif timepoint == '1F':
      ax.plot(dist_values, [x+5 for x in map_values],'--', label=timepoint, c='orange')
    elif timepoint == '1F_L98':
      ax.plot(dist_values, [x+5 for x in map_values], label=timepoint, c='orange')
    elif timepoint == '2F':
      ax.plot(dist_values, [x+8 for x in map_values],'--',label=timepoint, c='green', linewidth=linewidth)
    elif timepoint == '2F_L98':
      ax.plot(dist_values, [x+10 for x in map_values], label=timepoint, c='green', linewidth=linewidth)
    else:
      ax.plot(dist_values, map_values, label=timepoint)

plt.axvline(0.0, color='k', alpha=0.3)
plt.axvline(-4.8, color='k', alpha=0.3)
plt.axvline(-5.1, color='k', alpha=0.3)
plt.yticks([])
plt.xticks([])
plt.tight_layout()
plt.show()
#plt.savefig('MN1MN4_map_trace_1d_anom_fig2_lm.pdf')
