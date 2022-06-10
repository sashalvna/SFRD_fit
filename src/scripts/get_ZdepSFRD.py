import numpy as np
import astropy.units as u
import scipy


########################################################
## Chrip mass
########################################################
def Mchirp(m1, m2):
    chirp_mass = np.divide(np.power(np.multiply(m1, m2), 3./5.), np.power(np.add(m1, m2), 1./5.))
    return chirp_mass    


########################################################
## # SFR(z) Madau & Dickinson 2014 shape
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

    
########################################################
##  The mettalicity distribution dP/dZ(z)
########################################################
def skew_metallicity_distribution(max_redshift = 10.0,redshift_step = 0.01,
                                  mu_0=0.025, mu_z=-0.048, sigma_0=1.125, sigma_z=0.048,
                                  alpha = -1.767, 
                                  min_logZ  =-12.0, max_logZ  =0.0, step_logZ = 0.01,
                                  metals = [], redsh = [],
                                  min_logZ_COMPAS = np.log(1e-4),max_logZ_COMPAS = np.log(0.03)):
    #                                mu_0=0.025, mu_z=-0.048, sigma_0=1.125, sigma_z=0.048,alpha = -1.767,                               
    """
    Calculate the distribution of metallicities at different redshifts using a log skew normal distribution
    the log-normal distribution is a special case of this log skew normal distribution distribution, and is retrieved by setting 
    the skewness to zero (alpha = 0). 
    Based on the method in Neijssel+19. Default values of mu0=0.035, mu_z=-0.23, sigma_0=0.39, sigma_z=0.0, alpha =0.0, 
    retrieve the dP/dZ distribution used in Neijssel+19

    NOTE: This assumes that metallicities in COMPAS are drawn from a flat in log distribution!

    Args:
        max_redshift       --> [float]          max redshift for calculation
        redshift_step      --> [float]          step used in redshift calculation
        min_logZ_COMPAS    --> [float]          Minimum logZ value that COMPAS samples
        max_logZ_COMPAS    --> [float]          Maximum logZ value that COMPAS samples
        
        mu_0    = 0.025    --> [float]           location (mean in normal) at redshift 0
        mu_z    = -0.05    --> [float]           redshift scaling/evolution of the location
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
    
    ##################################
    # Follow Langer & Norman 2007? in assuming that mean metallicities evolve in z as:
    mean_metallicities = mu_0 * 10**(mu_z * redshifts) 
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
        step_logZ         = np.diff(log_metallicities)[0]
#         step_logZ         = step_logZ[0]
        print('step_logZ', step_logZ)
        
    ##################################
    # probabilities of log-skew-normal (without the factor of 1/Z since this is dp/dlogZ not dp/dZ)
    dPdlogZ = 2./(sigma[:,np.newaxis]) * normal_PDF((log_metallicities -  mu_metallicities[:,np.newaxis])/sigma[:,np.newaxis]) * normal_CDF(alpha * (log_metallicities -  mu_metallicities[:,np.newaxis])/sigma[:,np.newaxis] )

    ##################################
    # normalise the distribution over al metallicities
    norm = dPdlogZ.sum(axis=-1) * step_logZ #<< Fit does not converge if you multiply by step_logZ
    dPdlogZ = dPdlogZ /norm[:,np.newaxis]

    ##################################
    # assume a flat in log distribution in metallicity to find probability of drawing Z in COMPAS
    p_draw_metallicity = 1 / (max_logZ_COMPAS - min_logZ_COMPAS)
    
    return dPdlogZ, redshifts, metallicities, step_logZ, p_draw_metallicity


