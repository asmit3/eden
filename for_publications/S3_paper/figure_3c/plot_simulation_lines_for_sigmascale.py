from libtbx.easy_pickle import load
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
from scitbx.array_family import flex


# color scheme:
# https://matplotlib.org/stable/users/explain/colors/colormaps.html
mpl.rcParams['axes.linewidth'] = 2.0
data = load(sys.argv[1])
template = sys.argv[1].split('.pkl')[0]
all_r_mn, all_s00, all_s25, all_s50, all_s75, all_s100, dist_values, map_values = data 
timepoints = ['2F']
pair = sys.argv[2]
monomer='lm'
for timepoint in timepoints:
  map_values_sigma = flex.double()
  dist_values_sigma = flex.double()
  with open('sigma_scale_csv/%s_%s_%s_map_values_noNormalization.csv'%(timepoint,monomer, pair), 'r') as fin:
    for line in fin:
      if line !='\n':
        bx = line.split(',')
        if 'atoms' in line:
          atom1=bx[1].strip()
          atom2=bx[2].strip()
          continue
        if timepoint in ['2F_oldrestr', '2F']:
          map_values_sigma.append(float(bx[1]))
        else:
          map_values_sigma.append(float(bx[1]))
        dist_values_sigma.append(float(bx[0]))

z2=flex.linear_regression(map_values, map_values_sigma)
slope = z2.slope()
y_intercept = z2.y_intercept()
print ('Slope and y_intercept for sigma values vs absolute map values = %.5f, %.5f'%(slope, y_intercept))
scaled_s00 = [x*slope+y_intercept for x in all_s00]
scaled_s25 = [x*slope+y_intercept for x in all_s25]
scaled_s50 = [x*slope+y_intercept for x in all_s50]
scaled_s75 = [x*slope+y_intercept for x in all_s75]
scaled_s100 = [x*slope+y_intercept for x in all_s100]

# Actual plotting begins now
offset = 0 # 2.25 in large, 2.30 in sm
all_r_mn = [x-offset for x in all_r_mn]
dist_values = [x-offset for x in dist_values]
alpha=1.0
fig, ax = plt.subplots()
cmap=mpl.cm.get_cmap('pink')
cmap=cmap.reversed()
legend_text = 'Ox'
# linewidths
toy_lw = 2.5
data_lw = 7
ax.plot(all_r_mn, scaled_s00, label = f'%s=0%%'%legend_text, c=mpl.colors.rgb2hex(cmap(0.3), keep_alpha=True), linewidth=toy_lw, alpha=alpha)
ax.plot(all_r_mn, scaled_s25, label = f'%s=25%%'%legend_text, c=mpl.colors.rgb2hex(cmap(0.5), keep_alpha=True), linewidth=toy_lw, alpha=alpha)
ax.plot(all_r_mn, scaled_s50, label = f'%s=50%%'%legend_text, c=mpl.colors.rgb2hex(cmap(0.7), keep_alpha=True), linewidth=toy_lw, alpha=alpha)
ax.plot(all_r_mn, scaled_s75, label = f'%s=75%%'%legend_text, c=mpl.colors.rgb2hex(cmap(0.85), keep_alpha=True),linewidth=toy_lw, alpha=alpha)
ax.plot(all_r_mn, scaled_s100, label= f'%s=100%%'%legend_text, c=mpl.colors.rgb2hex(cmap(1.0), keep_alpha=True), linewidth=toy_lw, alpha=alpha)

axvline = {'data_lm_Mn1O6':1.75, 'data_lm_O5O6':2.1, 'data_lm_Ca1O6':2.70, 'data_sm_Mn1O6':1.75, 'data_sm_O5O6':2.0, 'data_sm_Ca1O6':2.60}
# Real data

ax.plot(dist_values, map_values_sigma,linewidth=data_lw, label='2F XFEL', c='forestgreen', alpha=alpha)
#ax.plot(dist_values_sigma, map_values_sigma, linewidth=2)
#ax2 = ax.twinx()
#ax2.plot(dist_values_sigma, map_values_sigma,linewidth=1, label='2F XFEL', c='forestgreen', alpha=alpha)
#plt.legend(fontsize=16, loc='upper left')
#plt.xlabel('Distances along O-Mn-O')
#plt.show()
plt.tight_layout()
plt.yticks([0, 5, 10, 15])
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
plt.axvline(axvline[template], alpha=0.3, c='grey')
plt.axvline(0.0, alpha=0.3, c='grey')
#plt.axvline(0.0, alpha=0.3, c='grey')
#plt.savefig('fig_%s.pdf'%template)
plt.show()



