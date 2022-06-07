
import h5py as h5
import numpy as np
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



############################
##PLOT setttings
from matplotlib import rc
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('font', family='serif')
rc('text', usetex=True)
fsize, SMALL_SIZE, MEDIUM_SIZE, BIGGER_SIZE = 30,20,20,30
for obj in ['axes','xtick','ytick']:
    plt.rc(obj, labelsize=MEDIUM_SIZE)          # controls default text sizes
for obj in ['figure','axes']:
    plt.rc(obj, titlesize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize



def Residuals_Zz_plane(obs_lookback = [], obs_metal = [], obs_SFRD = [],
                       model_t = None, model_redshift=None, model_y = None, model_SFRD = None,
                       chi_square_matrix = None, tmin = 0.0, tmax = 13.7, 
                       scatter_residuals = True, add_TNG = True, plot_dPdZcontours = True, 
                       neijssel_fit = True, dPdZ_text = None, SFR_text = None,
                       boundkleur = 'orange',COMPASkleur = 'Oranges', TNGkleur = 'YlGnBu', 
                       scatterKleur = 'RdYlGn_r', BBH_kleur = 'magma', savestr ='Residuals',
                       min_logZ_COMPAS = np.log(1e-4),max_logZ_COMPAS = np.log(0.03)
                      ):
    '''
    x, y, z             ---------------> redshift/lookback time, metallicities, dP/dZ
    tmin,tmax           ---------------> min and max time in Gyr to show as xlim 
    DCO, DCO_mask       ---------------> table of double compact objects, + mask of which to include in plot 
    kleur, kleurlabel   ---------------> colour/colour label of contour
    savestr             ---------------> string added to save name of plot
    min_logZ_COMPAS     ---------------> min ln(metal) that ocurs in COMPAS
    max_logZ_COMPAS     ---------------> max ln(metal) that ocurs in COMPAS
    '''
    ######################################
    # Create the Figure
    #     fig, ax = plt.subplots(figsize = (15,10))
    fig = plt.figure(figsize = (15,15))
    gs = gridspec.GridSpec(10, 8)
    gs.update(wspace=0.0, hspace=0.5)
    ax        = plt.subplot(gs[0:4, 1:7])
#     ax.set_facecolor('lightgrey')
    ax_metals = plt.subplot(gs[5:, :3])
    ax_redsh  = plt.subplot(gs[5:, 5:8])
    ######################################
    
    
    ##############################################################################
    # Load TNG data (either original or interpolated)
    if len(obs_SFRD) == 0:
        print('Using original TNG')
        with h5.File(TNGlocation+"SFRMetallicityFromGasTNG100.hdf5", "r") as f:
            MetalBins     = f["MetalBins"][:]
            Obs_Lookbacktimes = f["Lookbacktimes"][:]
            BoxSfr        = f["Sfr"][:]

        # Take the centers of the metallicity bins
        Obs_center_Zbin = (MetalBins[:-1] + MetalBins[1:])/2. 

        # Convert SFR from sfr/box to sfr Mpc-3
        littleh = 0.6774
        Rbox    = 75/littleh
        Obs_cosmic_SFR = BoxSfr / Rbox**3 *u.Mpc**-3
        Obs_cosmic_SFR = Obs_cosmic_SFR.value
        Obs_cosmic_SFR = Obs_cosmic_SFR.T

        ##########################################
        # "observed" TNG metallicities that we use for our calculations
        ##########################################
        center_Zbin = (MetalBins[:-1] + MetalBins[1:])/2.
        # Let's not use ALL metallicities in the TNG.. (they go waay too low!)
        low_bound_Z_ind = np.where(center_Zbin > 1e-5)[0]# index of center_Zbin, where Z > 1e-5
        # Let's not use ALL metallicities in the TNG.. (they go waay too low!)
        bound_Z_ind = np.where(np.logical_and(center_Zbin > 1e-5, center_Zbin < 50*0.014))[0]# index of center_Zbin, where Z < 50 Zsun
        tofit_TNG_metals = center_Zbin[bound_Z_ind]   


    else:
        print('Using interpolated TNG')
        Obs_Lookbacktimes = obs_lookback
        Obs_center_Zbin   = obs_metal
        Obs_cosmic_SFR    = obs_SFRD
        
    # Convert observed lookback times to observed redshifts
    Obs_redshifts     = [z_at_value(cosmo.lookback_time,t*u.Gyr) for t in obs_lookback[1:]] 
    Obs_redshifts.insert(0,0) # put redshift zero at the start
    
    ##############################################################################
    # Top panel: SFRD
    ##############################################################################
    ######################################
    # plot TNG
    if add_TNG:   
        ######################################
        # now actually plot it
        tng_color_map = sns.light_palette("#fe875d", as_cmap=True) # construct smooth cmap from one colour
        tng_colors   = tng_color_map(np.linspace(0.,1.0, 7) )      # Pick 7 colours from this smooth colourmap
        tng_color     = ListedColormap(tng_colors)                 # Turn it back into a cmap

        TNG = ax.pcolormesh(Obs_Lookbacktimes, Obs_center_Zbin, Obs_cosmic_SFR, 
                            rasterized=True, norm=matplotlib.colors.LogNorm(vmin=1e-8,vmax=1e-1), 
                            cmap=tng_color, alpha=0.95 ) #matplotlib.cm.YlGnBu
        cbaxes1 = fig.add_axes([0.925, 0.1, 0.03, 0.8]) #[left, bottom, width, height]
        cb = plt.colorbar(TNG, cax = cbaxes1, label= r"$\mathrm{TNG \ SFRD \ [M_{\odot} yr^{-1} Mpc^{-3}]}$")  


        
    ##############################################################################
    # YOUR dP/dZ MODEL
    ##############################################################################
    if plot_dPdZcontours:
        # Plot the contours of the metallicity density
        levels = [1e-7,1e-6, 1e-5, 1e-4,1e-3,1e-2,5e-1]#np.logspace(2, 10., 10+1)
        COMPAS_cmap = sns.color_palette(COMPASkleur, as_cmap=True)
        
        cs = ax.contour(model_t, model_y, model_SFRD, levels, linewidths=4, cmap=COMPAS_cmap,
                         locator=ticker.LogLocator(), alpha = 0.95, zorder=0)
        ax.clabel(cs,inline=1,fontsize=20, levels = levels, use_clabeltext=True, fmt = '%.0e')

        ###################
        # print all fit values
        ax.text(0.82, 0.88, '$\mathrm{Fit \ parameters:}$',  fontsize=22, transform=plt.gcf().transFigure)
        # dPdZ
        ax.text(0.82, 0.75, dPdZ_text, fontsize=20, transform=plt.gcf().transFigure,
               bbox=dict(facecolor='none', edgecolor='grey', boxstyle='round,pad=0.5'))
        # SFR
        ax.text(0.82, 0.63, SFR_text, fontsize=20, transform=plt.gcf().transFigure,
               bbox=dict(facecolor='none', edgecolor='grey', boxstyle='round,pad=0.5'))
    
    ###################
    #Plotvalues
    ax.xaxis.grid(5) # vertical lines
    ax.set_yscale('log')
    ax.set_xlabel('$\mathrm{Lookback \ time \ [Gyr]}$', fontsize = 25)
    ax.set_ylabel('$\mathrm{Metallicity}$', fontsize = 25)
    
    ######################################
    #### Add redshift Axis ####
    ax2 = ax.twiny()

    redshift_tick_list = [0,0.1, 0.25, 0.5, 0.75, 1.0,1.5, 2, 3, 6, 10]
    # Find loockback location for each of our ages
    z_ticks = [cosmo.lookback_time(z) for z in redshift_tick_list]
    
    # And annotate the tick labels :)
    ax2.set_xticks([cosmo.lookback_time(z).value for z in redshift_tick_list])
    ax2.set_xticklabels(['${:g}$'.format(z) for z in redshift_tick_list],Fontsize = 20)
    ax2.set_xlabel('$\mathrm{redshift}$', fontsize = 25)

    #Make sure top and bottom axis are lined up (have same limmits)
    ax.set_xlim(tmin, tmax)
    ax2.set_xlim(tmin, tmax)
    
    ax.set_ylim(1e-5, 1e0)
    
    ##############################################################################
    # COMPA fiducial
    ##############################################################################
    if neijssel_fit:
        # We don't have that many metallicities in TNG, which is ugly. 
        # to solve this I'm making a new list of metallicities that contains the TNG fit metallicities
        n_increase = 2
        log_tofitmetals = np.log10(tofit_TNG_metals) # TNG data is equally spaced in log
        dlog_fit_metal = np.diff(log_tofitmetals)[0]/n_increase # we divide the step size by n_increase
        neijssel_metallicities = 10**np.arange(min(log_tofitmetals), max(log_tofitmetals)+dlog_fit_metal, dlog_fit_metal )


        # Some check for closeness to original?
        print(np.diff(np.log10(tofit_TNG_metals)) )
        log_tofitmetals = np.log10(tofit_TNG_metals)
        dlog_fit_metal = np.diff(log_tofitmetals)[0]/2
        new_array = 10**np.arange(min(log_tofitmetals), max(log_tofitmetals)+dlog_fit_metal, dlog_fit_metal )
        print(tofit_TNG_metals, new_array)

        print(tofit_TNG_metals-new_array[::2] < 1e-10 )

        if (tofit_TNG_metals-new_array[::n_increase] > 1e-10 ).any():
            raise ValueError('your new metal array does not contain the old values')
        
        # Get dPdZ   #         neijssel_metals = np.logspace(-5., -0.5, 50)
        neijssel_dPdlogZ, neijssel_redshifts, neijssel_metallicities, step_logZ, p_draw_metallicity = \
                        Z_SFRD.skew_metallicity_distribution(mu0=0.035, muz=-0.23,
                                                      alpha_0 = 0, alpha_z = 0, 
                                                      sigma_0=0.39, sigma_z =0, 
                                                      min_logZ  =-12.0, max_logZ  =0.0, step_logZ = 0.01,
                                                      metals=neijssel_metallicities, redsh = model_redshift)
        #Convert redshift to lookback time
        t_lookback = cosmo.lookback_time(neijssel_redshifts)
        #####################################
        # Get the SFR Neijssel et al 2019:
        neijssel_sfr = Z_SFRD.Madau_Dickinson2014(neijssel_redshifts, a=0.01, b=2.77, c=2.9, d=4.7) # Msun year-1 Mpc-3 
        Neijssel_SFRDzZ = neijssel_sfr*neijssel_dPdlogZ.T
        Neijssel_SFRDzZ = Neijssel_SFRDzZ.value
        
        cs_N = ax.contour(t_lookback, neijssel_metallicities, Neijssel_SFRDzZ, levels, linewidths=4,linestyles =':', alpha = 0.95, zorder=0,
                          cmap=sns.color_palette('Greys', as_cmap=True),locator=ticker.LogLocator())
        
    ##############################################################################
    # LEFT BOTTOM: SFRDz with metals on x-axis
    ##############################################################################
    redsfift_indces = [0,5,10,20, 40, 60,120,160]#[0,8,15,32, 50, 100, 150, 170]#
    colors     = plt.cm.coolwarm(np.linspace(0.,1.0, len(redsfift_indces))) #3rd num is the number of colours
    LAB = 'TNG'
    plot_lines = []
    
    # Plot a set of redshifts with offset
    for i, redshift_i in enumerate(redsfift_indces):
        if i != 0:
            LAB = None
        shift_step = (len(redsfift_indces)-1)*0.01 - 0.01*i
        ######################################
        # Observed: TNG data
        # print("!! error", "np.log10(Obs_center_Zbin)", np.log10(Obs_center_Zbin), "np.shape(Obs_cosmic_SFR)", np.shape(Obs_cosmic_SFR))
        # ax_metals.plot(np.log10(Obs_center_Zbin), Obs_cosmic_SFR[:,redshift_i] + shift_step,
                             # lw = 9, c = 'darkgrey', label = LAB)  

        ######################################
        # Model: new SFRD
        l = ax_metals.plot(np.log10(model_y), model_SFRD[:,redshift_i] + shift_step,
                       lw = 4, ls = '--', c = colors[i], label = "$z=%s$"%(np.round(Obs_redshifts[redshift_i], 3)) )    
        plot_lines.append([l])
        
        ######################################
        # Model: OLD (Neijssel et al. 2019)       
        if neijssel_fit:
            #print('Redshift used in Neijssel: ', np.round(neijssel_redshifts[redshift_i],3))
            ax_metals.plot(np.log10(neijssel_metallicities), Neijssel_SFRDzZ[:,redshift_i] + shift_step,
                             lw = 3, c = 'slategrey', ls = ':', zorder=0, alpha = 0.5)         
        
    ######################################
    # Plot values
    lines = ax_metals.get_lines()
    legend1 = ax_metals.legend(lines[0:1], ['TNG','TNG'],
                     bbox_to_anchor=(0.0, 1.0), loc='lower left')
    ax_metals.add_artist(legend1)
    
    legend2 = ax_metals.legend(lines[2:3], ['Neijssel et al. 2019',' c'],
                     bbox_to_anchor=(0.4, 1.0), loc='lower left')
    ax_metals.add_artist(legend2)
    
    legend3 = ax_metals.legend(lines[1::3], ["$%s$"%(np.round(Obs_redshifts[redshift_i], 3)) for redshift_i in redsfift_indces],
                     ncol=1, bbox_to_anchor=(1.02, 1), loc='upper left', title='$\mathrm{redshift}$')
    
    ax_metals.set_xlabel('$\log_{10}(Z)$', size =25)
    ax_metals.set_ylabel('$\mathrm{SFRD [yr^{-1}\ Mpc^{-3}}$]', size =25)

    ax_metals.set_ylim(-0.005, 0.08)
    ax_metals.set_xlim(-4, -0.5)
    
    
    ##############################################################################
    # Right BOTTOM: SFRDz with redshift on x-axis
    ##############################################################################
    metal_indices = [0,8,12,15,17,20,21,24]
    colors     = plt.cm.PiYG(np.linspace(0.,1.0, len(metal_indices))) #3rd num is the number of colours
    LAB = 'TNG'
    plot_lines = []
    
    # Plot a set of redshifts with offset
    for j, metal_i in enumerate(metal_indices):
        if i != 0:
            LAB = None
        shift_step = (len(metal_indices)-1)*0.01 - 0.01*j 
        ######################################
        # Observed: TNG data
        # ax_redsh.plot(Obs_redshifts, Obs_cosmic_SFR[metal_i,:] + shift_step,
                             # lw = 9, c = 'darkgrey', label = LAB)  

        ######################################
        # Model: NEW SFRD
        l = ax_redsh.plot(model_redshift, model_SFRD[metal_i,:] + shift_step,
                       lw = 4, ls = '--', c = colors[j] )    
        plot_lines.append([l])
        
        ######################################
        # Model: OLD (Neijssel et al. 2019)       
        if neijssel_fit:
            #print('Metal used in Neijssel: ', np.log10(neijssel_metallicities[metal_i*n_increase]) )
            ax_redsh.plot(neijssel_redshifts, Neijssel_SFRDzZ[metal_i*n_increase,:] + shift_step,
                             lw = 3, c = 'slategrey', ls = ':', zorder=0, alpha = 0.5)      
        
    ######################################
    # Plot values
    lines = ax_redsh.get_lines()
    legend1 = ax_redsh.legend(lines[0:1], ['TNG','TNG'],
                     bbox_to_anchor=(0.0, 1.0), loc='lower left')
    ax_redsh.add_artist(legend1)
    
    legend2 = ax_redsh.legend(lines[2:3], ['Neijssel et al. 2019',' c'],
                     bbox_to_anchor=(0.4, 1.0), loc='lower left')
    ax_redsh.add_artist(legend2)
    
    legend3 = ax_redsh.legend(lines[1::3], ["$%s$"%(np.round(np.log10(Obs_center_Zbin[metal_i]), 1)) for metal_i in metal_indices],
                     ncol=1, bbox_to_anchor=(1.02, 1), loc='upper left', title='$\log_{10}Z$')
    
    ax_redsh.set_xlabel('$\mathrm{redshift}$', size =25)
    ax_redsh.set_ylabel('$\mathrm{SFRD [yr^{-1}\ Mpc^{-3}}$]', size =25)
    
    ax_redsh.set_ylim(-0.005, 0.075)
    
    
    ##############################################################################
    print('saving here', paths.figures / 'SFRD_FIT_evaluation_compare.pdf')
    fig.savefig(paths.figures / 'SFRD_FIT_evaluation_compare.pdf',  bbox_inches='tight', dpi=300)
    # plt.savefig(save_loc + 'SFRD_FIT_evaluation_compare_1.pdf',  bbox_inches='tight')
    
    # plt.show()




#####################################
#####################################
mu0_best     = 0.025
muz_best     = -0.049
sigma0_best  = 1.129
sigmaz_best  = 0.048
alpha0_best  = -1.778
alphaz_best  = 0.0


# lets interpolate at regular redshift intervals
z_new    = np.arange(0, 10.1, 0.05)                       # new redshifts
xnew     = [cosmo.lookback_time(z).value for z in z_new]  # corresp lookback times
ynew     = np.logspace(-5., -0.5, 100) # new metals

#####################################
#####################################

#####################################
# Get dPdZ 
# Make sure to calculate it at the same redshifts/metals as used in your fit
dPdlogZ, redshifts, metallicities, step_logZ, p_draw_metallicity = \
                Z_SFRD.skew_metallicity_distribution(mu0=mu0_best, muz=muz_best,alpha_0 = alpha0_best, alpha_z = alphaz_best, 
                                              sigma_0=sigma0_best, sigma_z =sigmaz_best, min_logZ  =-12.0, max_logZ  =0.0, step_logZ = 0.01,
                                              metals=ynew, redsh = z_new)
#Convert redshift to lookback time
t_lookback = cosmo.lookback_time(redshifts)


#####################################
#####################################
sf_a_best     = 0.017
sf_b_best     = 1.481
sf_c_best     = 4.452
sf_d_best     = 5.913

#####################################
#####################################

# Get the SFR
sfr = Z_SFRD.Madau_Dickinson2014(redshifts, a=sf_a_best, b=sf_b_best, c=sf_c_best,  d=sf_d_best) # Msun year-1 Mpc-3 
MSSFR = sfr*dPdlogZ.T
#####################################

fit_values_string = '$\mathrm{log-skew-normal}$'+'\n'+\
                    '$\mu_0=%s,$'%(np.round(mu0_best,3)) +'\n'+\
                    '$\mu_z=%s,$'%(np.round(muz_best,3)) +'\n'+\
                    '$\sigma_0=%s,$'%(np.round(sigma0_best,3)) +'\n'\
                    '$\sigma_z=%s,$'%(np.round(sigmaz_best,3)) +'\n'\
                    '$a_0=%s$'%(np.round(alpha0_best,3))

SFR_fit_string = '$\mathrm{Star \ formation \ rate}$'+'\n'+\
                 '$a=%s,$'%(np.round(sf_a_best,3)) +'\n'+\
                 '$b=%s,$'%(np.round(sf_b_best,3)) +'\n'+\
                 '$c=%s,$'%(np.round(sf_c_best,3)) +'\n'\
                 '$d=%s,$'%(np.round(sf_d_best,3))

Residuals_Zz_plane(obs_lookback = xnew, obs_metal = ynew, obs_SFRD = [],
                   model_t = t_lookback.value, model_redshift=redshifts, model_y = metallicities, model_SFRD =MSSFR.value,
                   #chi_square_matrix = chi_square_zZ,
                   scatter_residuals = False, add_TNG = False, 
                   plot_dPdZcontours = True, neijssel_fit = True,
                   COMPASkleur="crest", dPdZ_text = fit_values_string, SFR_text = SFR_fit_string)#light:k#404040 fe1100


