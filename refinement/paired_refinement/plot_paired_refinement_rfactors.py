import sys
from libtbx.easy_pickle import load
import matplotlib.pyplot as plt
import numpy as np

fname = sys.argv[1]
# Input should be the .pickle file written out by paired_refinement.py
data = load(fname)

delta_rwork = data['delta_rwork']
delta_rfree = data['delta_rfree']
error_rwork = data['error_rwork']
error_rfree = data['error_rfree']
delta_rwork_randomized=data['delta_rwork_randomized']
delta_rfree_randomized=data['delta_rfree_randomized']
resolution_bins=data['resolution_bins']


indices = np.arange(len(delta_rwork))
width=0.15
plt.bar(indices, delta_rwork, width, label='delta Rwork')
plt.bar(indices+width, delta_rfree, width, label='delta Rfree')
plt.bar(indices+2*width, delta_rwork_randomized, width, label='delta Rwork randomized', hatch='//')
plt.bar(indices+3*width, delta_rfree_randomized, width, label='delta Rfree randomized', hatch='//')
plt.bar(indices+4*width, [0.0], width)
plt.legend()
plt.ylabel('R-factor')
plt.xlabel('Resolution bin')
plt.xticks(indices, resolution_bins, rotation=45)

plt.show()

