# import necessary libraries to run choppe and analyze its results

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd

# set the plot style

params= { 'font.family':'sans-serif',
            "font.weight":"bold",
                'xtick.labelsize':10,
                'ytick.labelsize':10
    }

matplotlib.rcParams.update(params)

# read and store the input spectra
# these are spectra of Jupiter collected by Juno/JIRAM during its first perijove passage
# they are a fraction (40 %) of the ones published in Biagiotti et al. 2025 "Machine-Learning Spectral Clustering Techniques: Application to Jovian Clouds from Juno/JIRAM and JWST/NIRSpec"
# they are also available at: https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/JUNO/jiram.html
# if you wish to use them outside this example you are gently asked to contact the JIRAM PI Alessandro Mura at alessandro.mura@inaf.it

spectra_file=pd.read_csv("input_spectra_40.txt",  header=None, sep='\s+', engine='python')
spectra=spectra_file.to_numpy()

# import chopper and set the flags

import chopper as chopper

labels, clusters, coefficients, mean_spectra, mean_std= chopper.do_Clustering(input_spectra=spectra, flag_solar_correction=False, flag_PCA=False,
                                                      flag_GMM=False, flag_N_clusters=False)


# plot the results 

K=len(mean_spectra[:,0])
wv=np.linspace(2, 5, len(spectra[0]))


for i in range(K):
    plt.plot(wv, mean_spectra[i], label="Cluster N. "+str(int(i+1)))
    plt.fill_between(wv, mean_spectra[i]-mean_std[i], mean_spectra[i]+mean_std[i], alpha=0.4)
plt.xlim(min(wv), max(wv))
plt.xlabel(r"Wavelength ($\mu$m)", fontweight="bold")
plt.ylabel(r"Radiance (W sr$^{-1}$ m$^{-2}$ $\mu$m$^{-1}$)", fontweight="bold")
plt.legend()
plt.show()
#plt.savefig("images//mean_cluster_spectra.pdf", dpi=500)
#plt.close()

# plot the results using an ancillary file
# these are the longitude and latitude values associated to the Jupiter spectra collected by Juno/JIRAM during PJ1 
# they are published in Biagiotti et al. 2025 "Machine-Learning Spectral Clustering Techniques: Application to Jovian Clouds from Juno/JIRAM and JWST/NIRSpec"
# they are also available at: https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/JUNO/jiram.html
# if you wish to use them outside this example you are gently asked to contact the JIRAM PI Alessandro Mura at alessandro.mura@inaf.it

# in particular in this code we replicate and A.7

ancillary_file=pd.read_csv("input_log_40.txt",  header=None, sep='\s+', engine='python')

ancillary1=ancillary_file[0]
ancillary2=ancillary_file[1] 

label_ax="System III West Longitude (°)"
label_ay="Planetocentric Latitude (°)"

ca1=[]
ca2=[]

for i in range(K):
    ca1.append([])
    ca2.append([])

for i in range(K):
    ca1[i].append([ancillary1[j] for j in clusters[i]])
    ca2[i].append([ancillary2[j] for j in clusters[i]])

dimensions=len(coefficients[0])

clabels=[]

for i in range(dimensions):
    clabels.append("c"+str(int(i+1)))

num_variables=dimensions
num_plots = num_variables
fig, axes = plt.subplots(num_plots, num_plots, figsize=(10, 10))
fig.subplots_adjust(hspace=0.05, wspace=0.05)
for idx, ax in enumerate(axes.flat):
    i = idx // num_variables 
    j = idx % num_variables   

    if j > i:
        ax.set_visible(False)
        continue

    if i == j:
        for cas in range(K):
            x_data=[coefficients[bix,j] for bix in clusters[cas]]
            ax.hist(x_data, bins=50, density=True, histtype='stepfilled', alpha=0.7)
            ax.set_yticks([]) 
    else:
        for hera in range(K):
            x_data=[coefficients[kanan,j] for kanan in clusters[hera]]
            y_data=[coefficients[kanan,i] for kanan in clusters[hera]]
            ax.plot(x_data, y_data, '.', markersize=1, alpha=0.3)


    if i != num_variables - 1:
        ax.set_xticks([])
    else:
        ax.set_xlabel(clabels[j], fontsize=14, fontweight="bold")
    if j != 0:
        ax.set_yticks([])
    else:
        ax.set_ylabel(clabels[i], fontsize=14, fontweight="bold")
    if i == j and i != num_variables - 1:
        ax.set_xticks([])
    ax.tick_params(direction='in', top=True, right=True, labelsize=10)

plt.show()
#plt.savefig("images//cornerplot.pdf", dpi=500)
#plt.close()


for i in range(K):
    plt.scatter(ca1[i], ca2[i], s=80, label="Cluster N. "+str(int(i+1)))
plt.gca().invert_xaxis()
plt.xlabel(label_ax, fontweight="bold")
plt.ylabel(label_ay, fontweight="bold")
plt.legend()
plt.show()
#plt.savefig("images//ancillary_distribution.pdf", dpi=500)
#plt.close()