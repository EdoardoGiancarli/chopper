# chopper.py
![Chopper Logo](images/chopper_logo.png)
Hi fellow rebels :) Welcome to the chopper.py repository

## What is chopper.py
`chopper.py` is a code developed during my PhD in order to automatically find similar groups of spectra in a generic collection. I'm an astronomer therefore I used it mainly to study the Jupiter spectra acquired by the NASA Juno mission and by JWST during its Eearly Release Science phase. If you're interested you can find the details in a forthcoming paper currently under review in Astronomy and Astrophysics.

However, I'm not an imperial, therefore you can use for any spectral collection ! More in general also for any collection of N elements each of them described by M features.

The code is very simple and here there is a diagram showing what it does:
![Chopper flowchart](images/chopper_flowchart.png)

## How to install

```
git clone https://github.com/astro-francy/chopper.git
cd chopper
pip install requirements.txt
```

## How to use
I suggest you to write a run.py or run.ipynb file awhere you read or import your spectral collection `spectra`. If so then you need to just implement these lines:
```
import chopper as chopper

labels, clusters, coefficients, mean_spectra, mean_spectra_std= chopper.do_Clustering(input_spectra=spectra, flag_solar_correction=False, flag_PCA=False,
                                                      flag_GMM=False, flag_N_clusters=False)
```
Let's wrap up the inputs and outputs of this function. 

### `input_sectra`
This is your spectral collection. It should be a `np.ndarray` with shape `(N,M)`

### `flag_solar`
```
labels, clusters, coefficients, mean_spectra= chopper.do_Clustering(input_spectra=spectra, flag_solar_correction=False, flag_PCA=False,
                                                      flag_GMM=False, flag_N_clusters=False, au=5.2, incidence_angle=inc, wavelengths=wv)
```
This is the first `bool` flag to be checked. It is particularly useful if your spectral collection regards solar system objects. If `flag_solar=False` the code does not modify your spectra. If `flag_solar=True` it modifies your spectra in BRDF values assuming a Lambertian model (https://en.wikipedia.org/wiki/Lambert%27s_cosine_law). To do that you also need to specify three further input parameters:
1. au (float): the distance of the object from the sun in A.U.
2. incidence_angle (ndarray): a N-dimensional array containing incidence angles in (°)
3. wavelengths (ndarray); a M-dimensional array containing the wavelength values expressed in micron

N.B. chopper.py was designed to be a general tool to take a first look at large spectral collections. Therefore it approximates the solar irradiance as a black body. This is a good approximation in most of the cases.

### 'flag_PCA'
```
labels, clusters, coefficients, mean_spectra= chopper.do_Clustering(input_spectra=spectra, flag_solar_correction=False, flag_PCA=True,
                                                      flag_GMM=False, flag_N_clusters=False)
```
This is the second 'bool' flag to be checked. If 'flag_PCA=False' chopper retains only the PCA components that cumulatively explain more than 95 % of the dataset variance. If 'flag_PCA=True' chopper displays the first 'components_to_be_displayed' (line 63 of chopper.py) components and let the user decide through an input which variance threshold to consider.

### 'flag_GMM'
```
labels, clusters, coefficients, mean_spectra= chopper.do_Clustering(input_spectra=spectra, flag_solar_correction=False, flag_PCA=False,
                                                      flag_GMM=True, flag_N_clusters=False)
```
This is the third 'bool' flag to be checked. If 'flag_GMM=False' chopper performs GMM clustering on the PCA decomposition coefficients distribution obtained after checking the second flag or not. If 'flag_GMM=True' chopper displays the components retained at the previous step and the relative decomposition coefficients. Then the user decide through an input which coefficients to use or not to perform the GMM clustering.

### 'flag_N_clusters'
```
labels, clusters, coefficients, mean_spectra= chopper.do_Clustering(input_spectra=spectra, flag_solar_correction=False, flag_PCA=True,
                                                      flag_GMM=False, flag_N_clusters=True)
```
This is the last 'bool' flag to be checked. If 'flag_N_clusters=False' chopper uses the elbow method in order to choose the number of clusters 'K' to use during the GMM clustering. If 'flag_N_clusters=True' chopper let the user decide through an input which 'K' to consider. 

### outputs
1. 'labels': N-dimensional np.array associating any element in 'spectra' to a number between 1 and 'K'. i.e. the resulting cluster
2. 'clusters': list composed by 'K' lists. Each of them conatins the indexes (referred to 'spectra') of the elements composing each retrieved cluster
3. 'coefficients': np.ndarray containing the PCA decomposition coefficients
4. 'mean_spectra': np.ndarray containing the mean spectra computed for each cluster
5. 'mean_spectra_std': np.ndarray containing the standard deviation associated to each cluster mean spectrum

# Examples
In the directories 'grs', 'ganymede', and 'raman' there are some examples on how to build a run.py file and how to plot and analyze the chopper outputs. In particular these are:

1. grs: Contains Juno/JIRAM spectra of the Great Red Spot Region on Jupiter. The notebook replicates Fig. 9 and A7 from Biagiotti et al. 2025 (submitted to A&A). In particular it explains how to use both 'flag_PCA' and 'flag_GMM'
2. ganymede: Contains Juno/JIRAM spectra of Ganymede. The notebook shows how chopper.py is able to automatically replicate the results of https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020JE006508. It also explains how it can be useful to use 'flag_solar_correction'
3. raman: it requires ramanspy (https://ramanspy.readthedocs.io/en/latest/) to be installed. The notebook shows how chopper.py is able to automatically clusterize non-astrophysical spectra of bacteria. Is it also useful to show how 'flag_N_clusters' can be used.

# Notes

v1 : uploaded in June 2025. Basic code as described in Biagiotti et al. 2025 (submitted to A&A)
