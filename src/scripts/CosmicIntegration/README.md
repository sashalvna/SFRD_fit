
<h1>
Cosmic integration
</h1>
<p>

The scripts in this folder will perform the cosmic integration on $\texttt{COMPAS\\_Output\\_wWeights.h5}$ and produce a file called $\texttt{Rate\\_info.h5}$.
These files can be downloaded from <a href="https://zenodo.org/record/7612755">Zenodo.</a> 

If you rename the files (in particular $\texttt{COMPAS\\_Output\\_wWeights.h5}$ ), make sure to update <a href="../init_values.py"> init_values.py</a> correspondingly.

</p>
<h3> Step 1: run CallCosmicIntegration.py </h3>
<p>
The only file you have to worry about is <a href="./CallCosmicIntegration.py"> CallCosmicIntegration.py</a>, which is intended to run in a slurm based HPC environment.  Make sure to change the 'SlurmJobString' to whatever suits you environment, or at the very least change the 'partitions' flag to the name of whichever partition you would like to use. 

The flag `ZdepSFRD_param_sets' determines the shape of your cosmic starformation history, $\mathcal{S}(z,Z)$. It takes an array of arrays: [[dP/dZ], [SFRD(z)]] (the first containing 5 metallicity distributionparameters [mu0, muz, sigma0, sigmaz, alpha0], and the second the 4 SFR parameters [sf_a, sf_b, sf_c, sf_d]).

<ahref="./CallCosmicIntegration.py">CallCosmicIntegration.py</a> calls <a href="./FastCosmicIntegration.py"> FastCosmicIntegration.py</a>, which is where the real magic happens (see also <a href="https://github.com/TeamCOMPAS/COMPAS/blob/dev/online-docs/notebooks/CosmicIntegration.ipynb"> here</a> for an extended tutorial on how this works under the hood). 

To succesfully run <a href="./FastCosmicIntegration.py"> FastCosmicIntegration.py</a>, you need all of the following files to be present: ClassCOMPAS.py
SNR_Grid_IMRPhenomPv2_FD_all_noise.hdf5, selection_effects.py, and totalMassEvolvedPerZ.py.

On full resolution, every cosmic integration run job takes between 10 - 45 min, depending on your variation. 
 
</p>
<h3> Step 2: run CheckCompletionAndCombine.py </h3>

<p>
The previous step outputs a file called CI_job_IDs.txt. 
If you now run <a href="./CheckCompletionAndCombine.py"> CheckCompletionAndCombine.py</a>, it will collect all the job IDs from the cosmic integration, check for complettion, and combine all the individual hdf5 files into Rate_info.h5 .

All done :)!
</p>


