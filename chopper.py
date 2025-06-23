import sys
import corner
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture as GMM
from sklearn.cluster import KMeans
from kneed import KneeLocator


def do_Clustering(input_spectra, flag_solar_correction, flag_PCA, flag_GMM, flag_N_clusters, distance_in_au=1, incidence_angle=0, wavelengths=[2]):

    print("#########################################################################")
    print("################## directly from the ghost crew...... ###################")
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
    print("################# chopper.py v1 written by @astro-francy ################")
    print("#########################################################################")

    if not isinstance(input_spectra, np.ndarray):
        raise TypeError("input_spectra should be a ndarray !!!")
    
    if not np.issubdtype(input_spectra.dtype, float):
        raise TypeError("input_spectra is not composed entirely by numbers !!!")
    
    if not isinstance(flag_solar_correction, bool):
        raise TypeError("flag_solar_correction should be either True or False !!!")
    
    if not isinstance(flag_PCA, bool):
        raise TypeError("flag_PCA should be either True or False !!!")
    
    if not isinstance(flag_GMM, bool):
        raise TypeError("flag_GMM should be either True or False !!!")
    
    if not isinstance(flag_N_clusters, bool):
        raise TypeError("flag_N_clusters should be either True or False !!!")


    params= { 'font.family':'sans-serif',
            "font.weight":"bold",
                'xtick.labelsize':10,
                'ytick.labelsize':10
    }

    matplotlib.rcParams.update(params)
    components_to_be_displayed=10
    labelx="Spectral Index"
    labely="Spectral magnitude"


    spectra=input_spectra
    nspectra=len(spectra[:,0])

    features=len(spectra[0])  

    wv=np.linspace(1, features, features)

    if flag_solar_correction==True:
        #solar=pd.read_csv("solar_irradiance.txt",  header=None, sep='\s+', engine='python', skiprows=1)
        if not isinstance(incidence_angle, np.ndarray):
            raise TypeError("incidence_angle should be a ndarray !!!") 
        if len(incidence_angle)!=nspectra:
            raise TypeError("incidence_angle has not the right dimensions !!!") 
        if not isinstance(wavelengths, np.ndarray):
            raise TypeError("wavelengths should be a ndarray !!!") 
        if len(wavelengths)!=features:
            raise TypeError("wavelengths has not the right dimensions !!!") 
        if not isinstance(distance_in_au, float):
            raise TypeError("au should be a float !!!") 

        def solar_irradiance(wavelength_um, temperature_k):
            h = 6.62607015e-34  
            c = 2.99792458e8    
            k_B = 1.380649e-23  
            wavelength_m = wavelength_um * 1e-6

            exp_term = np.exp((h * c) / (wavelength_m * k_B * temperature_k))

            radiance = (2 * h * c**2) / (wavelength_m**5 * (exp_term - 1))

            return radiance* 1e-6 * (0.0046505**2) * np.pi

        au=distance_in_au
        inc=incidence_angle
        wv2=wavelengths
        factor_bdrf=solar_irradiance(wv2,5777)/(np.pi * (au**2))
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

    if flag_PCA==True:

        if not isinstance(labelx, str):
            raise TypeError("labels should be str !!!") 
        
        if not isinstance(labely, str):
            raise TypeError("labels should be str !!!") 

        print("You checked the PCA flag.")
        print("Now I will display you the first "+str(components_to_be_displayed)+" PCA eigenvectors and the relative cumulative variance explained by each of them")
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
        plt.close()

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


    if flag_GMM==True:
        print("You checked the #1 GMM flag.")
        print("Now I will display you the PCA eigenvectors")
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
        plt.close()

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



    if flag_N_clusters==True:
        print("You checked the #2 GMM flag.")

        question3=input("Please insert the desired number of cluster:")
        q3=str(question3)

        if q3=="":
            sys.exit("Error: You did not insert any threshold !!!")   
        elif int(q3)<2:
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

        '''
        plt.plot(range(2, 11), wcss)
        plt.axvline(x=suggested_n_clusters, color="k", ls="--")
        plt.title('Elbow Method')
        plt.xlabel('Number of clusters', fontweight="bold")
        plt.ylabel('WCSS metric', fontweight="bold")
        plt.savefig("images//kelbow.pdf", dpi=400)
        plt.cla()
        '''


        if suggested_n_clusters:
            print(f"Based on k-elbow method, a possible optimal number of clusters is: {suggested_n_clusters}")
        else:
            sys.exit("It's difficult to automatically determine a clear elbow point from the data. Please use the flag and inspect the plot.")

        K=suggested_n_clusters


    gmm = GMM(n_components=K, random_state=42)
    gmmlabels=gmm.fit(data).predict(data)

    clusters=[]
    cspectra=np.zeros((K, len(spectra[0])))
    cstd=np.zeros((K, len(spectra[0])))

    for i in range(K):
        clusters.append([])

    for i in range(len(gmmlabels)):
        for j in range(K):
            if gmmlabels[i]==j:
                clusters[j].append(i)


    for i in range(K):
        cspectra[i,:]=np.mean([spectra[j] for j in clusters[i]], axis=0)
        cstd[i,:]=np.std([spectra[j] for j in clusters[i]], axis=0)

    
    print("")
    print("/~/_|  |_\~\ --- Done. Bye :)")
    print("")

    return gmmlabels, clusters, data, cspectra, cstd