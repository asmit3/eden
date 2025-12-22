from __future__ import absolute_import, division, print_function
from cctbx import maptbx
from scipy.optimize import curve_fit
import numpy as np

from scitbx.array_family import flex
import matplotlib as mpl
    
    
def gauss(x, amplitude, mu, sigma):
  return amplitude*np.exp(-(x-mu)*(x-mu)/(2*sigma*sigma))

def triple_gauss(x, amplitude1, mu1, sigma1, amplitude2, mu2, sigma2, amplitude3, mu3, sigma3):
  return amplitude1*np.exp(-(x-mu1)*(x-mu1)/(2*sigma1*sigma1)) + amplitude2*np.exp(-(x-mu2)*(x-mu2)/(2*sigma2*sigma2)) + amplitude3*np.exp(-(x-mu3)*(x-mu3)/(2*sigma3*sigma3)) 

# https://stackoverflow.com/questions/25668828/how-to-create-colour-gradient-in-python
def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

def closest_index(arr, value):
  """
  Returns the index of the element in `arr` that is closest to `value`,
  or None if `value` is outside the [min(arr), max(arr)] range.

  No assumptions about arr.
  """
  if not arr: return None
  min_val = min(arr)
  max_val = max(arr)
  if value < min_val or value > max_val: return None
  closest_i = 0
  min_diff = abs(arr[0] - value)
  for i in range(1, len(arr)):
    diff = abs(arr[i] - value)
    if diff < min_diff:
      min_diff = diff
      closest_i = i
  return closest_i

def run():
  all_r_mn = []
  all_vals_mn = []
  all_r_o = []
  all_vals_o = []
  all_r_o2 = []
  all_vals_o2 = []
  all_s = []
  #

  all_s00 = []
  all_s25 = []
  all_s50 = []
  all_s75 = []
  all_s100 = []

  # Mn 43, q0=1.00 works quite well

  monomer='lm'
  b_iso_mn = {'lm':44, 'sm':46} # 45 for sm
  b_iso_o = {'lm':28, 'sm':30} # 24.5 for sm
  b_iso_o2 = {'lm':30, 'sm':30}
  d_min = 1.93
  radius = 5


  mn = maptbx.atom_curves(scattering_type="Mn", scattering_table="wk1995")
  o  = maptbx.atom_curves(scattering_type="O", scattering_table="wk1995")
  im_mn  = mn.image(d_min=d_min, b_iso=b_iso_mn[monomer], radius_min=-radius, radius_max=radius)
  im_o   =  o.image(d_min=d_min, b_iso=b_iso_o[monomer], radius_min=-radius, radius_max=radius)
  im_o2  =  o.image(d_min=d_min, b_iso=b_iso_o2[monomer], radius_min=-radius, radius_max=radius)
  baseline_distance1 = {'lm':2.25, 'sm':2.30}
  baseline_distance2 = {'lm':4.35, 'sm':4.25}
  dist1  = baseline_distance1[monomer] # LM  = 2.25
  dist2  = baseline_distance2[monomer] # LM  = 4.35
  q0 = 1.00
  q1 = 1.00
  q2 = 0.65
  for r_mn, vals_mn,  r_o, vals_o,  r_o2, vals_o2  in zip(im_mn.radii,  im_mn.image_values,
                                                          im_o.radii,   im_o.image_values,
                                                          im_o2.radii,  im_o2.image_values):

    i = closest_index(arr=im_o.radii+dist1, value=r_mn)

    s = q0*vals_mn
    s00 = q0*vals_mn
    s25 = q0*vals_mn
    s50 = q0*vals_mn
    s75 = q0*vals_mn
    s100 = q0*vals_mn
    if i is not None:
      s = s + q1*im_o.image_values[i]
      s00 = s00 + q1*im_o.image_values[i]
      s25 = s25 + q1*im_o.image_values[i]
      s50 = s50 + q1*im_o.image_values[i]
      s75 = s75 + q1*im_o.image_values[i]
      s100 = s100 + q1*im_o.image_values[i]

    i = closest_index(arr=im_o2.radii+dist2, value=r_mn)
    if i is not None:
      s = s + q2*im_o2.image_values[i]
      s00 = s00 + 0.00*im_o2.image_values[i]
      s25 = s25 + 0.25*im_o2.image_values[i]
      s50 = s50 + 0.50*im_o2.image_values[i]
      s75 = s75 + 0.75*im_o2.image_values[i]
      s100 = s100 + 1.00*im_o2.image_values[i]

    all_r_mn.append(r_mn)
    all_vals_mn.append(vals_mn)
    all_r_o.append(r_o+dist1)
    all_vals_o.append(vals_o)
    all_r_o2.append(r_o2+dist2)
    all_vals_o2.append(vals_o2)
    all_s.append(s)
    # Add all the population sums here
    all_s00.append(s00)
    all_s25.append(s25)
    all_s50.append(s50)
    all_s75.append(s75)
    all_s100.append(s100)

  # Fit gaussians here
  #p0= [4.0, 0.00, 0.80, 1.5, 1.75, 0.80, 1.5, -2.0, 0.80]  
  #params, pcov= curve_fit(triple_gauss, all_r_mn, all_s, p0=p0)
  map_values = flex.double() 
  dist_values = flex.double()
  offset = {'lm':2.1, 'sm':2.15}
  with open('absolute_scale_csv/2F_%s_O5O6_map_values_noNormalization.csv'%monomer, 'r') as fin:
    for line in fin:
      if line !='\n':
        bx = line.split(',')
        if 'atoms' in line:
          atom1=bx[1].strip()
          atom2=bx[2].strip()
          continue
        map_values.append(float(bx[1]))
        dist_values.append(float(bx[0])+offset[monomer]) # Subtracting here to center the plot

  import matplotlib.pyplot as plt
  fig, ax = plt.subplots()
  #ax.scatter(all_r_mn, all_vals_mn, c='red', label='Mn only')
  #ax.scatter(all_r_o, all_vals_o,c='royalblue', label='O only')
  #ax.scatter(all_r_o2, all_vals_o2, c='silver',label='O2 only')
  #ax.scatter(all_r_mn, all_s, label = 'Sum of map values')


  c1='#1f77b4' #blue
  c2='green' #green


  cmap=mpl.cm.get_cmap('Purples')


  ax.scatter(all_r_mn, all_s00, label = 'Ox=0.00', c=mpl.colors.rgb2hex(cmap(0.3), keep_alpha=True), linewidth=1)
  ax.scatter(all_r_mn, all_s25, label = 'Ox=0.25', c=mpl.colors.rgb2hex(cmap(0.5), keep_alpha=True), linewidth=1)
  ax.scatter(all_r_mn, all_s50, label = 'Ox=0.50', c=mpl.colors.rgb2hex(cmap(0.7), keep_alpha=True), linewidth=1)
  ax.scatter(all_r_mn, all_s75, label = 'Ox=0.75', c=mpl.colors.rgb2hex(cmap(0.85), keep_alpha=True), linewidth=1)
  ax.scatter(all_r_mn, all_s100, label= 'Ox=1.00', c=mpl.colors.rgb2hex(cmap(1.0), keep_alpha=True), linewidth=1)

  #ax.scatter(all_r_mn, all_s00, label = 'Sum of map values, Ox=0.00', c=colorFader(c1, c2, 0.1), linewidth=1)
  #ax.scatter(all_r_mn, all_s25, label = 'Sum of map values, Ox=0.25', c=colorFader(c1, c2, 0.25), linewidth=1)
  #ax.scatter(all_r_mn, all_s50, label = 'Sum of map values, Ox=0.50', c=colorFader(c1, c2, 0.50), linewidth=1)
  #ax.scatter(all_r_mn, all_s75, label = 'Sum of map values, Ox=0.75', c=colorFader(c1, c2, 0.75), linewidth=1)
  #ax.scatter(all_r_mn, all_s100, label = 'Sum of map values, Ox=1.00', c=colorFader(c1, c2, 1.0), linewidth=1)
  # Real data
  ax.scatter(dist_values, map_values, label='2F XFEL data', c='orange')


  # Now show the fits
  if False:
    dist_values = all_r_mn
    ax.plot(all_r_mn, triple_gauss(all_r_mn, *params),'-o', c='gold', alpha=0.2)
    ax.plot(dist_values, gauss(dist_values, params[0], params[1], params[2]), '--', c='red', label='fit_Mn')
    ax.plot(dist_values, gauss(dist_values, params[3], params[4], params[5]), '--', c='royalblue', label='fit_O')
    ax.plot(dist_values, gauss(dist_values, params[6], params[7], params[8]), '--', c='silver', label='fit_O2')
    ax.annotate('%.2f, %.2f, %.2f'%(params[0], params[1], params[2]), xy=(params[1], params[0]), xytext=(params[1], params[0]))
    ax.annotate('%.2f, %.2f, %.2f'%(params[3], params[4], params[5]), xy=(params[4], params[3]), xytext=(params[4], params[3]))
    ax.annotate('%.2f, %.2f, %.2f'%(params[6], params[7], params[8]), xy=(params[7], params[6]), xytext=(params[7], params[6]))
  
  plt.axvline(dist1, alpha=0.2)
  plt.axvline(dist2, alpha=0.2)
  plt.legend()

  from libtbx.easy_pickle import dump
  dumper = [all_r_mn, all_s00, all_s25, all_s50, all_s75, all_s100, dist_values, map_values]
  dump('data_%s_O5O6.pkl'%monomer, dumper)

  plt.legend()
  #plt.xlabel('Distances along O-Mn-O')
  plt.tight_layout()
  plt.title('%.2f, %.2f, %.2f'%(b_iso_mn[monomer], b_iso_o[monomer], b_iso_o2[monomer]))
  #plt.savefig('fig_O5O6.pdf', format='pdf')
  plt.show()
  #from IPython import embed; embed(); exit()
if (__name__ == "__main__"):
  run()
