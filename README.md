
<h1>
An analytical fit to the metallicity-dependent cosmic star formation history, $\mathcal{S}(z,Z)$
</h1>
<p>
This is the source code associated to the paper titled:
</p>

<h3>
<a href="https://ui.adsabs.harvard.edu/abs/2022arXiv220903385V/abstract">The locations of features in the mass distribution of merging binary black holes are robust against uncertainties in the metallicity-dependent cosmic star formation history.</a>
</h3>

<p>

The files in src/scripts/ form the core of this work (i.e. reproduce all the Figures). 

Additionally, we provide a Jupyter notebook called <a href="./src/scripts/Notebooks/Fit_model_to_sfrdzZ.ipynb">Fit_model_to_sfrdzZ.ipynb</a> that fits our $\mathcal{S}(z,Z)$ model to an arbitrary input star formation rate - metallicity grid. Our fiducial is fit to data from the TNG 100 simulation ( $\texttt{SFRMetallicityFromGasTNG100.hdf5}$ ). 

  
The COMPAS binary population simulation data that was used to calculate the BBH mass distribution can be found <a href="https://sandbox.zenodo.org/deposit/1153294">on Zenodo.</a> If you would like to re-run the Figures in this work with your own variation of the cosmic starformation history, you can use the scripts in the <a href="./src/scripts/CosmicIntegration/">CosmicIntegration!</a>

If anything is unlcear, don't hesitate to shoot me a message, find my contact info  <a href="https://liekevanson.github.io/contact.html"> here</a>
  
</p>


<br>
<br>
<a href="https://github.com/LiekeVanSon/SFRD_fit/actions/workflows/build.yml">
<img src="https://github.com/LiekeVanSon/SFRD_fit/actions/workflows/build.yml/badge.svg?branch=main" alt="Article status"/>
</a>
<a href="https://github.com/LiekeVanSon/SFRD_fit/raw/main-pdf/arxiv.tar.gz">
<img src="https://img.shields.io/badge/article-tarball-blue.svg?style=flat" alt="Article tarball"/>
</a>
<a href="https://github.com/LiekeVanSon/SFRD_fit/raw/main-pdf/ms.pdf">
<img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
</a>
</p> 
This work was made using
<a href="https://github.com/showyourwork/showyourwork">
<img width = "150" src="https://raw.githubusercontent.com/showyourwork/.github/main/images/showyourwork.png" alt="showyourwork"/>
</a>
