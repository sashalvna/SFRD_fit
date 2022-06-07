import h5py as h5
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import astropy.units as u
from astropy import constants as const
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

from scipy import interpolate
from scipy.optimize import minimize
from scipy.optimize import curve_fit

from astropy.cosmology import WMAP9, z_at_value
from astropy.cosmology import Planck18  as cosmo# Planck 2018
from astropy.cosmology import z_at_value

# from astropy.table import Table, Column
import os

import paths

import matplotlib
from pylab import *
from matplotlib import ticker, cm
import matplotlib.gridspec as gridspec

def Mchirp(m1, m2):
    chirp_mass = np.divide(np.power(np.multiply(m1, m2), 3./5.), np.power(np.add(m1, m2), 1./5.))
    return chirp_mass    
   
base_dir    = '/Users/lieke/surfdrive/Documents/RateMassRedshiftEvolution/'
save_loc    =  base_dir+'/plots/'
TNGlocation = '/Users/lieke/surfdrive/Documents/CompareCOMPAS/'

############################
# Custom scripts
import get_ZdepSFRD as Z_SFRD



# SFR(z) Madau & Dickinson 2014 shape


########################################################
##
########################################################
def Madau_Dickinson2014(z, a=0.015, b=2.77, c=2.9, d=5.6):
    """
    Args:
        z             --> [list of floats] List of redshifts at which to calculate things
        a,b,c,d       --> [floats] values to determine the shape of our SFR
    
    Calculates the star-formation rate density as a function of redshift
    Based on the functional form from Madau & Dickinson 2014
    default 'Neijssel et al 2019': a=0.01, b=2.77, c=2.9,  d=4.7
    Madau & Dickinson 2014: a=0.015, b=2.7, c=2.9,  d=5.6
    Madau & Fragos 2017: a=0.01, b=2.6, c=3.2,  d=6.2

    Returns:
        SFR(z) in Msun/yr/Mpc^3
    """
    dm_dtdMpc = a * (1 + z)**b/( 1 + ( (1+z)/c )**d ) *u.Msun *u.yr**-1 *u.Mpc**-3
    return dm_dtdMpc # Msun year-1 Mpc-3 

    

# # Get dP/dZ(z)

# ### New log skewed distribution

# metallicity density distribution
# Basically this tells us what the probability of finding a certain metallicity at a certain redshift is. 

# See also eq. 17 in : https://arxiv.org/ftp/arxiv/papers/1501/1501.02344.pdf


# # (log) Skew normal dist: 

# PDF of the skew normal is defined as:

# \begin{equation}
# f(Z) = 2 \phi \left(\frac{Z - \mu}{\sigma}\right) \Phi\left(\alpha \frac{Z - \mu}{\sigma} \right)
# \end{equation}

# $\phi(t)$ and $\Phi(t)$ are the standard ($\sigma = 1$) normal PDF and CDF respectively:

# The skew normal is a generalization of the normal distribution, that gets back to a normal distribution for $\alpha = 0$.


# ***
# The PDF for the **log-skew-normal** are then given by just substituting the random variable Z, with ln(Z)


# \begin{equation}
# f(Z) = 2 \phi(\frac{ln(Z) - \mu}{\sigma}) \Phi(\alpha \frac{ln(Z) - \mu}{\sigma})
# \end{equation}

# or written out explicitely:

# \begin{equation}
# f(Z) = \frac{2}{Z \sigma \sqrt{2 \pi}} e^{\frac{-1}{2} \left(\frac{\ln(Z) - \mu}{\sigma}\right)^2} 
# \int_{-\infty}^{x = \left(\alpha \frac{\ln(Z) - \mu}{\sigma} \right) } \frac{1}{\sqrt{2 \pi}} e^{\frac{-1}{2}t^2} dt
# \end{equation}

# or

# \begin{equation}
# f(Z) = \frac{2}{Z \sigma \sqrt{2 \pi}} e^{\frac{-1}{2} \left(\frac{\ln(Z) - \mu}{\sigma}\right)^2} 
# \frac{1}{2} \left[ 1 + erf\left(\frac{ \left(\alpha \frac{\ln(Z) - \mu}{\sigma} \right)  }{\sqrt{2}} \right) \right]
# \end{equation}


# The extra $1/Z$ factors in front of the PDFs come from 
# \begin{equation}
# \frac{dP}{dZ} = \frac{dP}{d\ln Z} \frac{d\ln Z}{dZ} = \frac{dP}{d\ln Z} \frac{1}{Z}
# \end{equation}

# ***

# ## Moments
# Moments of a probability distribution: the zeroth moment is the total probability (i.e. one), the first moment is the expected value, the second central moment is the variance, the third standardized moment is the skewness, and the fourth standardized moment is the kurtosis. 

# equation 23) and 24) from : https://arxiv.org/ftp/arxiv/papers/1501/1501.02344.pdf
# give the mean and variance of our log-skew-normally distributed random variable Z:


# \begin{equation}
# E(Z) = 2 e^{\mu} e^{\sigma^2 /2} \Phi(\beta \sigma)
# \end{equation}


# \begin{equation}
# VAR(Z) = 2 e^{2\mu} e^{\sigma^2} (e^{\sigma^1} \Phi(2 \beta \sigma)  - 2\Phi(\beta \sigma)^2)
# \end{equation}


# with 
# \begin{equation}
# \beta = \frac{\alpha}{\sqrt{1 + \alpha^2} }
# \end{equation}



def skew_metallicity_distribution(max_redshift = 10.0,redshift_step = 0.01,
                                  mu0=0.025, muz=-0.048, sigma_0=1.125, sigma_z=0.048,
                                  alpha_0 = -1.767, alpha_z = 0, 
                                  min_logZ  =-12.0, max_logZ  =0.0, step_logZ = 0.01,
                                  metals = [], redsh = [],
                                  min_logZ_COMPAS = np.log(1e-4),max_logZ_COMPAS = np.log(0.03)):
    #                                mu_0=0.025, muz=-0.048, sigma_0=1.125, sigma_z=0.048,alpha = -1.767,                               
    """
    Calculate the distribution of metallicities at different redshifts using a log skew normal distribution
    the log-normal distribution is a special case of this log skew normal distribution distribution, and is retrieved by setting 
    the skewness to zero (alpha = 0). 
    Based on the method in Neijssel+19. Default values of mu0=0.035, muz=-0.23, sigma_0=0.39, sigma_z=0.0, alpha =0.0, 
    retrieve the dP/dZ distribution used in Neijssel+19

    NOTE: This assumes that metallicities in COMPAS are drawn from a flat in log distribution!

    Args:
        max_redshift       --> [float]          max redshift for calculation
        redshift_step      --> [float]          step used in redshift calculation
        min_logZ_COMPAS    --> [float]          Minimum logZ value that COMPAS samples
        max_logZ_COMPAS    --> [float]          Maximum logZ value that COMPAS samples
        
        mu0    = 0.025    --> [float]           location (mean in normal) at redshift 0
        muz    = -0.05    --> [float]           redshift scaling/evolution of the location
        sigma_0 = 1.25     --> [float]          Scale (variance in normal) at redshift 0
        sigma_z = 0.05     --> [float]          redshift scaling of the scale (variance in normal)
        alpha   = -1.77    --> [float]          shape (skewness, alpha = 0 retrieves normal dist)

        min_logZ           --> [float]          Minimum logZ at which to calculate dPdlogZ (influences normalization)
        max_logZ           --> [float]          Maximum logZ at which to calculate dPdlogZ (influences normalization)
        step_logZ          --> [float]          Size of logZ steps to take in finding a Z range

    Returns:
        dPdlogZ            --> [2D float array] Probability of getting a particular logZ at a certain redshift
        metallicities      --> [list of floats] Metallicities at which dPdlogZ is evaluated
        p_draw_metallicity --> float            Probability of drawing a certain metallicity in COMPAS (float because assuming uniform)
    """
    import scipy

    ##################################
    # the PDF of a standard normal distrtibution
    def normal_PDF(x):
        return 1./(np.sqrt(2* np.pi)) * np.exp(-(1./2) * (x)**2 )

    ##################################
    # the CDF of a standard normal distrtibution
    def normal_CDF(x):
        return 1./2. * (1 + scipy.special.erf(x/np.sqrt(2)) )
    
    ##################################
    if len(redsh) == 0:
        # Make redshifts
        redshifts = np.arange(0, max_redshift + redshift_step, redshift_step)
    else:
        redshifts = redsh
        
    ##################################
    # Experiment with redshift dependence sigma
    # LOG-LINEAR
    sigma = sigma_0*10**(sigma_z*redshifts)
    #  LINEAR   sigma = sigma_z*redshifts + sigma_0
    
        ##################################
    # Experiment with redshift dependent alpha
    # LINEAR
#     alpha = alpha_0 + (alpha_z*redshifts)
    # LOG-LINEAR (better)
    alpha = alpha_0*10**(alpha_z*redshifts)
    
    ##################################
    # Follow Langer & Norman 2007? in assuming that mean metallicities evolve in z as:
    mean_metallicities = mu0 * 10**(muz * redshifts) 
    #print('np.shape(mean_metallicities)', np.shape(mean_metallicities))
        
    # Now we re-write the expected value of ou log-skew-normal to retrieve mu
    beta = alpha/(np.sqrt(1 + (alpha)**2))
    PHI  = normal_CDF(beta * sigma) # phi is now sigma x alpha dimentional 
    mu_metallicities = np.log(mean_metallicities/(2.*PHI) * 1./(np.exp(0.5*sigma**2) )  ) 

    ##################################
    if len(metals) == 0:
        # create a range of metallicities (thex-values, or raandom variables)
        log_metallicities = np.arange(min_logZ, max_logZ + step_logZ, step_logZ)
        metallicities = np.exp(log_metallicities)
    else: 
        #use a pre-determined array of metals
        metallicities     = metals
        log_metallicities = np.log(metallicities)
        step_logZ         = np.diff(log_metallicities)
        step_logZ         = step_logZ[0]
        #print('step_logZ', step_logZ)
        
    ##################################
    # probabilities of log-skew-normal (without the factor of 1/Z since this is dp/dlogZ not dp/dZ)
    dPdlogZ = 2./(sigma[:,np.newaxis]) * normal_PDF((log_metallicities -  mu_metallicities[:,np.newaxis])/sigma[:,np.newaxis]) * normal_CDF(alpha[:,np.newaxis] * (log_metallicities -  mu_metallicities[:,np.newaxis])/sigma[:,np.newaxis] )

    ##################################
    # normalise the distribution over al metallicities
    norm = dPdlogZ.sum(axis=-1) #* step_logZ << Fit does not converge if you multiply by step_logZ
    dPdlogZ = dPdlogZ /norm[:,np.newaxis]

    ##################################
    # assume a flat in log distribution in metallicity to find probability of drawing Z in COMPAS
    p_draw_metallicity = 1 / (max_logZ_COMPAS - min_logZ_COMPAS)
    
    return dPdlogZ, redshifts, metallicities, step_logZ, p_draw_metallicity


