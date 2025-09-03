import numpy as np
import matplotlib.pyplot as plt
from IPython import embed
from iotbx import mtz
from scitbx.array_family import flex
import pandas as pd
import sys


seed = np.random.randint(0, 2**32 - 1)  
print("Using seed: ", seed)
np.random.seed(seed)


try: 
    user_input = str(input("Enter a seed (or press Enter to use random seed): ")).strip()
except SyntaxError as e: 
    print("No input detected. Using random seed:", seed)

# Initialize RNG with that seed

mtz_obj = mtz.object("/Users/yyklab/Desktop/eden/model_analysis_scripts/error_analysis/0F.mtz")
# /Users/kanishkkondaka/Desktop/phenix/0F.mtz 


dataset = mtz_obj.as_miller_arrays_dict()




key = ('crystal', 'Experimental-data-used-in-refinement', 'F-obs-filtered') # Fmodel 
#miller_array = dataset[key]
fobs_array = dataset[key]


key2 = ('crystal', 'Model-structure-factors-(bulk-solvent-and-all-scales-included)', 'F-model')
array = dataset[key2]
fmodel_array = dataset[key2]
fmodel_values = np.array(array.amplitudes().data())

data_values = fobs_array.data()

f_values = np.array(list(data_values))


sigmas = fobs_array.sigmas()
sigf_values = np.array(list(sigmas))

difference_f_values = np.abs((f_values - fmodel_values))/ np.abs(f_values)
difference = np.abs((f_values - fmodel_values))
delta_f = np.abs(f_values - fmodel_values)
overall_rwork = np.sum(np.abs((f_values - fmodel_values)))/(np.sum(f_values))
print("The overall Rwork is: ", overall_rwork)
R_work = np.sum(np.abs(f_values - fmodel_values))/np.sum(f_values)


#embed()

"""
1. Undersatnd the fobs code, i don't understand none
2.Protocol 2: Check that we are sumtarcting F_o(hkl) - F_c(hkl). we can check Fo.get.matching_indices, first try this, 
then do the other thing.  we can check that the miller index is indetical using for loop. make sure that you are plotting Fo - Fc, 
be careful. We need to calculate Fo - F/Fo which is called Rwork, which is also called reliability factors
2. when plotting resolution with prot 2, the change should be larger as the resolutioion increases. 
3. we can get resolution from code that asmit will give me that comes from F0, unit cell
4. Also calculate Rwork which is the sum of (Fo - Fc)/F over all indices. 

Other thigns: 
Friedels law is when Fo(h) = Fo(-h) 
this is broken when therei s metalloenzymes which have anomolos scatters
Fo(h) = fo + F' + if''
fo normal scattering, anomalous f' if''
"""

"""
Checking that hkl indices are the same
for i in range(len(f_values)): 
    if fobs_array.indices()[i] != fmodel_array.indices()[i]: 
        print("FALSE")

print("TRUE")
"""

uc = fobs_array.unit_cell()

d = uc.d(fobs_array.indices())

resolution_values = np.array(list(d))

fobsf = fobs_array.resolution_filter(d_min=2.08, d_max=2.312).data()

fobsf = np.array(list(fobsf))

fmodelf = fmodel_array.resolution_filter(d_min = 2.08, d_max = 2.312).amplitudes().data() # you need the .amplitudes
# .data() alone will just return complex numbers

fmodelf = np.array(list(fmodelf))

Rworkf = np.sum(np.abs(fobsf - fmodelf))/np.sum(fobsf)
print("The Rwork for last bin is: ", Rworkf)

"""
#Plotting Fobs vs Fmodel

slope, intercept = np.polyfit(f_values, fmodel_values, deg=1)
print("slope, intercept", slope, intercept)
plt.figure(figsize =(10,4))
plt.scatter(f_values, fmodel_values)
plt.title("Fc vs Fo")
plt.xlabel("Fc")
plt.ylabel("Fo")
plt.grid(True)
plt.show()
"""



#Plotting Resolution vs deltaf_over_fo
plt.figure(figsize =(10,4))
plt.scatter(resolution_values, difference/f_values)
#plt.ylim(0, 1)
plt.axhline(1.0)
plt.title('Resolution vs Rwork ')
plt.xlabel("Resolution")
plt.ylabel("Rwork")
plt.grid(True)
plt.gca().invert_xaxis()




df = pd.DataFrame({
    "resolution": resolution_values,
    "fobs": f_values,
    "fmodel": fmodel_values,
    "diff": difference
})

num_bins = 30
df['res_bin'], bin_edges = pd.qcut(df['resolution'], q=num_bins, retbins=True, labels=False, duplicates='drop')
#Quantile-based discretization function.
#Discretize variable into equal-sized buckets based quantiles
"""
bin_edges  [ 2.08000004  2.10467417  2.13006903  2.15654932  2.18413277  2.21321401
  2.24379646  2.27625119  2.31056444  2.347009    2.38571566  2.42747752
  2.47208242  2.5201006   2.57206988  2.62837804  2.69008766  2.75792086
  2.83328787  2.91715925  3.01228758  3.12098199  3.24690682  3.39638399
  3.57674719  3.80357712  4.10044342  4.51921889  5.1842151   6.55114551
 33.63764466]
"""

binned_means = np.array([])
weights = np.array([])
binned_numerators = [] # these arrays are for verification of Rwork value
binned_denominators = []
for i in range(num_bins):
    bin_df = df[df['res_bin'] == i]
    # selects all rows from the DataFrame df where the value in the 'res_bin' column is equal to i.
    # res_bin has all of the bins, so this extracts all the values in that bin and puts it into bin_df
    
    numerator = bin_df['diff'].sum()
    denominator = bin_df['fobs'].abs().sum()
    binned_denominators.append(denominator)
    binned_numerators.append(numerator)
    weight = len(bin_df) / len(df)
    # the weight is calculated as a proportion of the numbers that are in the bin over the length of the entire
    # dataframe
    weights = np.append(weights, weight)
    mean_value = (numerator / denominator) 
#    binned_means.append(mean_value)
    binned_means = np.append(binned_means, mean_value)
# bin width = difference between bin edges
bin_widths = np.diff(bin_edges)
# the widths are calculated from the difference between each bin. Since there is a large gap between 33 and 6, 
# it ends up being the biggest bin. This means that entire bin contributes the same amount to the total Rwork value
# as just one bin.
bin_starts = bin_edges[:-1]


rwork_direct = df['diff'].sum() / df['fobs'].abs().sum()

rwork_from_bins = np.sum(binned_numerators) / np.sum(binned_denominators)

print("Rwork (direct):", rwork_direct)
print("Rwork (from bins):", rwork_from_bins)

plt.figure(figsize=(10, 4))
plt.bar(bin_starts, binned_means, width=bin_widths, align='edge',
        color='skyblue', edgecolor='black')
plt.xlabel('Resolution')
plt.ylabel('Mean Rwork')
plt.title('Binned Rwork with Proportional Resolution Bins')
plt.gca().invert_xaxis()
plt.tight_layout()
plt.show()







r = np.random.choice(100, size=3, replace=False)
x_values = np.random.choice(len(f_values), size=3, replace=False)

# 
deltafobs_simulations = np.random.normal(loc=f_values[None, :], scale=delta_f[None, :], size=(100, len(f_values)))


plt.figure(figsize = (10,4))
for i, idx in enumerate(x_values, start=1):
    deltavalues = deltafobs_simulations[:, idx]
    deltavalues[deltavalues < 0] = 0
    plt.hist(deltavalues, bins=100, density=True, alpha=0.6, edgecolor='black', label='Fobs index {}'.format(idx))
plt.title("Protocol 2 PDF")
plt.ylabel("Probability Density")
plt.xlabel("Simulated Fobs")
plt.grid(True)


plt.figure(figsize=(10,4))
plt.xlabel("Fobs value")
plt.ylabel("Cumulative Probability")
#plt.title("CDFs of 3 of the Same Indices of 3 Different Simulations")
plt.title("Protocol 2 CDF")

for i, idx in enumerate(x_values, start=1):
    sorted_data = np.sort(deltafobs_simulations[:, idx])
    sorted_data[sorted_data < 0] = 0
    cdf = np.linspace(0, 1, len(sorted_data))
    plt.plot(sorted_data, cdf, label='Simulation {}'.format(i))








#outputting the 100 mtz files 
fobs_simulations = np.random.normal(loc=f_values[None, :], scale=sigf_values[None, :], size=(100, len(f_values)))

for i in range(len(fobs_simulations)): 
    x_flex = flex.double(fobs_simulations[i,:])

    new_miller_array = fobs_array.customized_copy(data=x_flex)

    mtz_dataset = new_miller_array.as_mtz_dataset(column_root_label="Fobs_perturbed")
    mtz_dataset.mtz_object().write("output"+ str(i+1)+".mtz")





sigma_over_f = pd.DataFrame({
    "resolution": resolution_values,
    "fobs": f_values,
    "sigf": sigf_values,
    "sigmaoverf": sigf_values/f_values
})

num_bins = 30
sigma_over_f['res_bin'], bin_edges = pd.qcut(sigma_over_f['resolution'], q=num_bins, retbins=True, labels=False, duplicates='drop')

binned_means = np.array([])
weights = np.array([])
binned_numerators = []
binned_denominators = []
for i in range(num_bins):
    bin_df = sigma_over_f[sigma_over_f['res_bin'] == i]

    binned_denominators.append(denominator)
    binned_numerators.append(numerator)
    weight = len(bin_df) / len(sigma_over_f)
    weights = np.append(weights, weight)
    mean_value = np.mean(bin_df["sigmaoverf"])
    binned_means = np.append(binned_means, mean_value)

bin_widths = np.diff(bin_edges)
bin_starts = bin_edges[:-1]

plt.figure(figsize=(10, 4))
plt.bar(bin_starts, binned_means, width=bin_widths, align='edge',
        color='skyblue', edgecolor='black')
plt.xlabel('Resolution ')
plt.ylabel('Fobs +/- SigF')
plt.title('Binned Resolution vs Fobs +/- SigF')
plt.gca().invert_xaxis()  
plt.tight_layout()
plt.show()





# plotting



plt.figure(figsize=(10, 4))
for i, idx in enumerate(x_values, start=1):
    values = fobs_simulations[:, idx]
    values[values < 0] = 0
    plt.hist(values, bins=100, density=True, alpha=0.6, edgecolor='black', label='Fobs index {}'.format(idx))

#plt.title("Value Distributions at 3 Fixed Fobs Indices of 3 Different Simulations")
plt.title("Protocol 1 PDF")
plt.xlabel("Simulated Fobs Value")
plt.ylabel("Probability Density")
plt.legend()
plt.grid(True)
plt.figure(figsize=(10, 4))
for i, idx in enumerate(x_values, start=1):
    sorted_data = np.sort(fobs_simulations[:, idx])
    sorted_data[sorted_data < 0] = 0
    cdf = np.linspace(0, 1, len(sorted_data))
    plt.plot(sorted_data, cdf, label='Simulation {}'.format(i))


plt.xlabel("Fobs value")
plt.ylabel("Cumulative Probability")
plt.title("Protocol 1 CDF")

plt.legend()
plt.grid(True)
plt.show()





"""
Rworkmodified.py 06/18/25

import numpy as np
import matplotlib.pyplot as plt
from iotbx import mtz
from scitbx.array_family import flex
import sys

fname = sys.argv[1]
mtz_obj = mtz.object(fname)
dataset = mtz_obj.as_miller_arrays_dict()


key = ('crystal', 'Experimental-data-used-in-refinement', 'F-obs-filtered') # Fmodel 
#miller_array = dataset[key]
fobs_array = dataset[key]

uc = fobs_array.unit_cell()
d = uc.d(fobs_array.indices())

key2 = ('crystal', 'Model-structure-factors-(bulk-solvent-and-all-scales-included)', 'F-model')
array = dataset[key2]
fmodel_array = dataset[key2]
fmodel_values = np.array(array.amplitudes().data())

data_values = fobs_array.data()

f_values = np.array(list(data_values))


sigmas = fobs_array.sigmas()
sigf_values = np.array(list(sigmas))
d_min=2.08 #2.08
d_max=2.312 #2.312


overall_rwork = np.sum(np.abs((f_values - fmodel_values)))/(np.sum(f_values))


delta_f_over_f = np.abs(f_values - fmodel_values)/f_values

import matplotlib.pyplot as plt
plt.figure()
plt.gca().invert_xaxis()
plt.scatter(d, delta_f_over_f)
plt.xlabel('Resolution(A)')
plt.ylabel('delta_f_over_f')
plt.axhline(1.0, c='k')

plt.figure()
plt.hist(delta_f_over_f, bins=100)
plt.xlabel('delta_f_over_f')
plt.show()


print ('Overall Rwork = %.2f'%overall_rwork)
fobsf = fobs_array.resolution_filter(d_min=d_min, d_max=d_max).data()
fmodelf = fmodel_array.resolution_filter(d_min = d_min, d_max = d_max).amplitudes().data()
bin_rwork = flex.sum(flex.abs(fobsf - fmodelf))/flex.sum(fobsf)
print('R-work value within resolution range %.3f and %.3f = %.2f'%(d_max, d_min, bin_rwork))


"""





