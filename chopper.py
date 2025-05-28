print("#########################################################################")
print("################## dierectly from the ghost crew...... ##################")
print("#########################################################################")

print("")
print("ˇ       --~~--_      ˇ ")
print("|____/~/_|  |_\~\____|" )
print("    |____________|                    Welcome user :)")
print("    |[][][][][][]|:=  .               I'm here to help you ")
print("  __| __         |__ \  ' .          / with your collection of spectra !")
print(" |  ||. |   ==   |  |  \    ' .     /  ") 
print("(|  ||__|   ==   |  |)   \      '<")
print(" |  |[] []  ==   |  |      \    '\|")
print(" |  |____________|  |        \    |")
print(" /__\     |_|    /__\          \ / \\")
print("")

print("#########################################################################")
print("################# chopper.py v0 written by @astro-francy ################")
print("#########################################################################")

print("")

print("")
print("/~/_|  |_\~\ --- Importing necessary libraries and setting plot style")
print("")

import sys
import corner
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture as GMM
from sklearn.cluster import KMeans
from kneed import KneeLocator
import csv
from scipy.interpolate import PchipInterpolator as CubicSpline
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize
import matplotlib.gridspec as gridspec
import matplotlib

params= { 'font.family':'sans-serif',
           "font.weight":"bold",
             'xtick.labelsize':10,
             'ytick.labelsize':10
}

matplotlib.rcParams.update(params)

components_to_be_displayed=10

print("")
print("/~/_|  |_\~\ --- Reading the files")
print("")

with open("chopper.txt", encoding='utf-8') as file:
    lines = file.readlines()

cl = [line.strip() for line in lines]

flag0=cl[3]
au=eval(cl[4])
flag1=cl[5]
flag2=cl[6]
flag3=cl[7]

labelx=cl[11]
labely=cl[12]

df=pd.read_csv(cl[0],  header=None, sep='\s+', engine='python')
spectra=df.to_numpy()
nspectra=len(spectra[:,0])

if cl[1]!="None":
    ancillary_number1=float(cl[14])
    ancillary_number2=float(cl[15])
    label_ax=cl[16]
    label_ay=cl[17]


    df2=pd.read_csv(cl[1], sep='\s+', header=None, engine='python')
    nspectra2=len(df2.index)

    if nspectra!=nspectra2:
        sys.exit("Error: Check your spectra or log pls ! They are not in the right dimensions")

    if ancillary_number1>max(df2.columns) or ancillary_number2>max(df2.columns):
        sys.exit("Error: Check the chopper.txt input file ! You are selecting a non existing column in your log file") 

    ancillary1=df2[ancillary_number1]
    ancillary2=df2[ancillary_number2]   

wv=np.linspace(float(cl[9]), float(cl[10]), len(spectra[0]))

if flag0=="y":
    solar=pd.read_csv("solar_irradiance.txt",  header=None, sep='\s+', engine='python', skiprows=1)
    if cl[1]=="None":
        sys.exit("Error: you didn't provide a log file, I cannot retrieve the incidence and emission angles !")  
    inc=df2[int(au[1])]
    csn = CubicSpline(solar[0].to_numpy()/1000, solar[1].to_numpy() * 1000)
    factor_bdrf=csn(wv)/(np.pi * (au[0]**2))
    i_on_f=[]
    for i in range(nspectra):
        angle=(inc[i]*np.pi)/180
        i_on_f.append(spectra[i]/(factor_bdrf*np.cos(angle)))
    spectra=np.array(i_on_f)

M=np.size(spectra[0])
mean=np.zeros(M)


for e in range(M):
   dummy=[]
   for i in range(nspectra):
      dummy.append(spectra[i][e])
   dummy=np.array(dummy)
   mean[e]=np.mean(dummy)


print("")
print("/~/_|  |_\~\ --- PCA step")
print("")

number_of_components=min(nspectra-1, len(spectra[0])-1 )

pca = PCA(int(number_of_components), random_state=0, svd_solver='full')
pca.fit(spectra)

evals = pca.explained_variance_ratio_
evals_cs = evals.cumsum()

evecs= np.vstack([pca.components_])

if flag1=="y":
    print("You checked the PCA flag.")
    print("Now I will display you the first"+str(components_to_be_displayed)+" PCA eigenvectors and the relative cumulative variance explained by each of them")
    print("Keep them in mind ^.^")

    fig, axs = plt.subplots(2,1,figsize=(8,7))

    for i in range(components_to_be_displayed):
        axs[0].plot(wv, evecs[i], label="PC"+str(int(i+1)))
    axs[0].set_xlabel(labelx, fontweight="bold")
    axs[0].set_ylabel(labely, fontweight="bold")
    axs[0].legend()
    axs[1].semilogx(np.arange(1, components_to_be_displayed+1, 1), evals_cs[:components_to_be_displayed], 'k.')
    axs[1].semilogx(np.arange(1, components_to_be_displayed+1, 1), evals_cs[:components_to_be_displayed], 'k-', alpha=0.7)
    axs[1].set_xlabel('PC Number', fontweight="bold")
    axs[1].set_ylabel('Cumulative Eigenvalues', fontweight="bold")

    plt.show()
    plt.cla()

    question1=input("Please insert the cumulative-variance threshold you desire to set:")
    q1=str(question1)

    if q1=="":
        sys.exit("Error: You did not insert any threshold !!!") 
    elif float(q1)>np.max(evals_cs):
        sys.exit("Error: the selected threshold is not reached !!!") 


    thr=float(q1)

else:
    thr=0.95

idx_thr=int(np.argmax(evals_cs>thr)+1)

if idx_thr==1:
    new_components=evecs[:2]
    dimensions=len(new_components)
else:
    new_components=evecs[:idx_thr]
    dimensions=len(new_components)

coeff=[]

for i in range(nspectra):
   coeff.append(np.dot(new_components, spectra[i] - mean))

c=np.array(coeff)

clabels=[]

for i in range(dimensions):
    clabels.append("c"+str(int(i+1)))

print("")
print("/~/_|  |_\~\ --- Clustering step")
print("")


if flag2=="y":
    print("You checked the #1 GMM flag.")
    print("Now I will display you the PCA eigenvectors again")
    print("Keep them in mind ^.^")

    for i in range(dimensions):
        plt.plot(wv, new_components[i], label="PC"+str(int(i+1)))
    plt.xlabel(labelx, fontweight="bold")
    plt.ylabel(labely, fontweight="bold")
    plt.legend()
    plt.show()  
    plt.close()  

    print("")
    print("Now I will display the corresponding PCA decomposition coefficients cornerplot")
    print("Keep it in mind ^.^") 

    fig=corner.corner(c, labels=clabels)
    plt.show()
    plt.cla()

    print("")

    print("Please insert the coefficients on which you desire that the GMM will be performed :)")
    question2=input("Please insert them in the form of a list (e.g. [3, 5, 6]) :")
    q2=np.array(eval(question2))-1
    dimensions=len(q2)

    if str(q2)=="":
        sys.exit("Error: You did not insert any number !!!") 
    elif len(q2)<2:
        sys.exit("Error: you selected only a single coefficient !!!") 

    dummy_arrays=[c[:,i] for i in q2]

    data=np.vstack(tuple(dummy_arrays)).T
    clabels=[clabels[i] for i in q2]

else:
    data=c



if flag3=="y":
    print("You checked the #2 GMM flag.")

    question3=input("Please insert the desired number of cluster:")
    q3=str(question3)

    if q3=="":
        sys.exit("Error: You did not insert any threshold !!!")   
    elif len(q3)<2:
        sys.exit("Error: you selected a number of clusters smaller than the number of the clustering dimensions !!!")  

    K=int(q3)

else:

    wcss=[]

    for i in range(2, 11):  # Try clustering from 1 to 10 clusters
        kmeans = KMeans(n_clusters=i, init='k-means++', max_iter=300, n_init=10, random_state=0, algorithm='elkan')
        kmeans.fit(data)
        wcss.append(kmeans.inertia_)

    kn = KneeLocator(range(2, 11), wcss, curve='convex', direction='decreasing')
    suggested_n_clusters=kn.knee

    plt.plot(range(2, 11), wcss)
    plt.axvline(x=suggested_n_clusters, color="k", ls="--")
    plt.title('Elbow Method')
    plt.xlabel('Number of clusters', fontweight="bold")
    plt.ylabel('WCSS metric', fontweight="bold")
    plt.savefig("images//kelbow.pdf", dpi=400)
    plt.cla()


    if suggested_n_clusters:
        print(f"Based on k-elbow method, a possible optimal number of clusters is: {suggested_n_clusters}")
    else:
        sys.exit("It's difficult to automatically determine a clear elbow point from the data. Please use the flag and inspect the plot.")

    K=suggested_n_clusters


gmm = GMM(n_components=K, random_state=42)
gmmlabels=gmm.fit(data).predict(data)


print("")
print("/~/_|  |_\~\ --- Plotting and printing the results")
print("")

clusters=[]
cspectra=np.zeros((K, len(spectra[0])))
ca1=[]
ca2=[]

for i in range(K):
    clusters.append([])
    ca1.append([])
    ca2.append([])

for i in range(len(gmmlabels)):
    for j in range(K):
        if gmmlabels[i]==j:
            clusters[j].append(i)


for i in range(K):
    cspectra[i,:]=np.mean([spectra[j] for j in clusters[i]], axis=0)
    ca1[i].append([ancillary1[j] for j in clusters[i]])
    ca2[i].append([ancillary2[j] for j in clusters[i]])


for i in range(K):
    plt.plot(wv, cspectra[i], label="Cluster N. "+str(int(i+1)))
plt.xlim(min(wv), max(wv))
plt.xlabel(labelx, fontweight="bold")
plt.ylabel(labely, fontweight="bold")
plt.legend()
plt.savefig("images//mean_cluster_spectra.pdf", dpi=500)
plt.close()

with open('outputs//mean_cluster_spectra.txt', 'w', newline="") as x:
    csv.writer(x, delimiter=" ").writerows(cspectra)


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
            x_data=[data[bix,j] for bix in clusters[cas]]
            ax.hist(x_data, bins=50, density=True, histtype='stepfilled', alpha=0.7)
            ax.set_yticks([]) 
    else:
        for hera in range(K):
            x_data=[data[kanan,j] for kanan in clusters[hera]]
            y_data=[data[kanan,i] for kanan in clusters[hera]]
            ax.plot(x_data, y_data, '.', markersize=1, alpha=0.3)


    if i != num_variables - 1:
        ax.set_xticks([])
    else:
        # Set x-labels for the bottom row
        ax.set_xlabel(clabels[j], fontsize=14, fontweight="bold")

    # Remove y-axis tick labels for all but the leftmost column
    if j != 0:
        ax.set_yticks([])
    else:
        # Set y-labels for the leftmost column
        ax.set_ylabel(clabels[i], fontsize=14, fontweight="bold")
    
    # Hide the x-axis ticks for the diagonal plots if they are not on the bottom row
    if i == j and i != num_variables - 1:
        ax.set_xticks([])

    # Fine-tune tick parameters for a cleaner look
    ax.tick_params(direction='in', top=True, right=True, labelsize=10)


plt.savefig("images//cornerplot.pdf", dpi=500)
plt.close()


if cl[1]!="None":

    log=[]

    for i in range(len(ancillary1)):
        log.append([ancillary1[i], ancillary2[i], int(gmmlabels[i]+1)])

    with open('outputs//clusters.txt', 'w', newline="") as x:
        csv.writer(x, delimiter=" ").writerows(log)

    for i in range(K):
        plt.scatter(ca1[i], ca2[i], label="Cluster N. "+str(int(i+1)))
    plt.xlabel(label_ax, fontweight="bold")
    plt.ylabel(label_ay, fontweight="bold")
    plt.legend()
    plt.savefig("images//ancillary_distribution.pdf", dpi=500)
    plt.close()

else:

    log=[]

    for i in range(len(ancillary1)):
        log.append([int(gmmlabels[i]+1)])

    with open('outputs//clusters.txt', 'w', newline="") as x:
        csv.writer(x, delimiter=" ").writerows(log)    


