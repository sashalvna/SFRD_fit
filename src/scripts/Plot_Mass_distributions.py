"""

# Plotting the BH mass distribution for several SFRD Z-distribution variations

"""

import sys
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
# configure backend here
matplotlib.use('Agg')

import seaborn as sns
from scipy import stats

import h5py as h5 
from astropy.table import Table
import astropy.units as u
from astropy import constants as const

# Chosen cosmology 
from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import z_at_value

import json
import argparse

# My own helper funcitons:
import importlib
import MassDistHelperFunctions as mfunc
importlib.reload(mfunc)

import gc
import paths


######################################
## locations
save_loc    =  str(paths.figures) + '/'
data_dir    =  str(paths.data) + '/'

rate_file       = '/RateData/Rate_info.h5'#'/RateData/small_Rate_info.h5' #
simulation_data = '/COMPAS_Output_wWeights.h5'#'/small_COMPAS_Output_wWeights.h5' #

only_stable = True 
only_CE = False
 

if np.logical_and(only_stable, only_CE):
    channel_string = 'all'
elif only_stable:
    channel_string = 'stable'
elif only_CE:
    channel_string = 'CE'

######################################
## PLOT setttings
plt.rc('font', family='serif')
from matplotlib import rc
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
fsize, SMALL_SIZE, MEDIUM_SIZE, BIGGER_SIZE = 30,25,25,30
for obj in ['axes','xtick','ytick']:
	plt.rc(obj, labelsize=MEDIUM_SIZE)          # controls default text sizes
for obj in ['figure','axes']:
	plt.rc(obj, titlesize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize




######################################
# Distribution plot function
######################################
def plot_mass_distribution(sim_dir = '', x_key = 'M_moreMassive', rate_keys = ['Rates_mu00.025_muz-0.05_alpha-1.77_sigma0%s_sigmaz0.05_zBinned'%(x) for x in [0.8, 1.125, 1.4]],
                   bins = np.arange(0.,55,2.5), z_bin_edges = [0,0.25], 
                   plot_LIGO = False, show_hist = False, show_KDE = True, kde_width = 0.1,  
                   only_stable = True, only_CE = True, 
                   bootstrap = False, bootstraps = 10, x_lim=(0.,50),  y_lim = (1e-2,30), 
                   Color = '#e388b0', linestyles = ['--','-', ':'], titletext = '',
                   labels = [r'$\mathrm{CE \ channel = \ }$', r'$\mathrm{stable \ RLOF \ channel = \ }$', r'$\mathrm{All = \ }$'],
                   xlabel = r'$M_{\mathrm{BH, 1}} \ \rm [M_{\odot}]$', ylabel = r'$\frac{d\mathcal{R}}{dM_{\mathrm{BH, 1} }} \ \mathrm{[Gpc^{-3}yr^{-1}M^{-1}_{\odot}]}$',
                   leg_args = {'loc':'lower left',  'bbox_to_anchor':[0.1, 0.], 'fontsize':20, 'title':''}, leg1_args = {'loc':'upper left', 'fontsize':18},
                   save_plot=False, save_name = 'Fiducial.png', multipanel = False, subplot = None):
    """
    Read DCO, SYS and merger rate data, necesarry to make the plots in this 

    Args:
        sim_dir              --> [string] Location of data

    Returns:
     plot

    """

    #########################################
    mass_binw = np.diff(bins)[0]

    plot_lines = []
    leg_labels = []

    #########################################
    # Start plotting
    if not multipanel: #(otherwise you define the Fig outside of this fucntion)
        fig, ax = plt.subplots(figsize = (12, 10))
    else:
        ax = subplot

    ################################################
    # GWTC-3 Powerlaw + Peak Mass distribution
    ################################################ 
    if plot_LIGO:
        color_plpeak = 'grey'#'#1f78b4'
        #################################################
        ## grab Powerlaw + Peak data from O3
        #################################################  
        #'/Volumes/StorageSpac/CompasOutput/output/'#'/n/home04/lvanson/LowMBH_peak/output/GWTC-3-population-data/analyses/PowerLawPeak/'
        input_fname = data_dir+'o3only_mass_c_iid_mag_iid_tilt_powerlaw_redshift_mass_data.h5'
        mass_1 = np.linspace(2, 100, 1000)
        mass_ratio = np.linspace(0.1, 1, 500)
        with h5.File(input_fname, "r") as f:
            mass_ppd = f["ppd"]
            mass_lines = f["lines"]
            mass_1_ppd = np.trapz(mass_ppd, mass_ratio, axis=0)
            mass_1_lower = np.percentile(mass_lines["mass_1"], 5, axis=0)
            mass_1_upper = np.percentile(mass_lines["mass_1"], 95, axis=0)
        ##############################
        # plot the max posterior and the 95th percentile
        ax.plot(mass_1, mass_1_ppd, lw=1.8, color=color_plpeak, zorder=1, label="$\mathrm{GWTC-3}$")
        ax.fill_between(mass_1, mass_1_lower, mass_1_upper, alpha=0.14,color=color_plpeak,zorder=0)

        legend1 = plt.legend(**leg1_args)

    nplot = 0

    ################################################
    # My Simulations
    ################################################
    try:
        DCO = mfunc.read_data(loc = sim_dir + '/' + simulation_data)
    except:
        print('data not found')

    # We'll show the change in rate between SFRD at several reference masses
    m10, m25, m40 = [], [], []

    print('nplot', nplot, '\n')
    ####################################################
    ### Loop over SFRD
    for i, rate_key in enumerate(rate_keys):
        print('rate_key', rate_key)

        # ### ## Reading Rate data ##
#         try:
            #     DCO_mask, redshifts, intrinsic_rate_density, intrinsic_rate_density_z0  = mfunc.read_rate_data(loc = sim_dir + '/Rate_info.hdf5', rate_key = rate_key)
        with h5.File(sim_dir + '/' + rate_file ,'r') as File:
            redshifts                 = File[rate_key]['redshifts'][()]
            # Different per rate key:
            DCO_mask                  = File[rate_key]['DCOmask'][()] # Mask from DCO to merging systems  
            #(contains filter for RLOF>CE and optimistic CE)
            intrinsic_rate_density    = File[rate_key]['merger_rate'][()]


#         except:
#             print('\n error reading', rate_key)
#             continue
            
        # # # # # # # # # # # # # # # # # # 
        #first bring it to the same shape as the rate table
        merging_BBH    = DCO[DCO_mask]

        #apply the additional mask based on your prefs
        if np.logical_and(only_stable, only_CE):
            print("Both only_stable and only_CE, I assume you just want both")
            channel_bool = np.full(len(merging_BBH), True)
        elif only_stable:
            channel_bool = merging_BBH['CE_Event_Count'] == 0
        elif only_CE:
            channel_bool = merging_BBH['CE_Event_Count'] > 0
        else:
            raise ValueError("Both only_stable =%s and only_CE=%s, set at least one to true"%(only_stable,only_CE))

        # we exclude CHE systems
        not_CHE = merging_BBH['Stellar_Type@ZAMS(1)'] != 16

        merging_BBH         = merging_BBH[not_CHE  * channel_bool]
        Red_intr_rate_dens  = intrinsic_rate_density[not_CHE * channel_bool, :]


        # # # # # # # # # # # # # # # # # # 
        ## Calculate average rate density per z-bin
        # crude_rate_density  = mfunc.get_crude_rate_density(Red_intr_rate_dens[:,:], redshifts, z_bin_edges)
        #########################################
        # X value and weight
        x_vals              = merging_BBH[x_key]
        i_redshift = np.where(redshifts == 0.2)[0][0]
        print('i_redshift', i_redshift)
        Weights             = Red_intr_rate_dens[:, i_redshift]#crude_rate_density[:,0]
        print('!X!X!X!', np.shape(Weights), np.shape(x_vals) )
        print(labels[i], ' len(table)=', len(merging_BBH) , ' Rate = ', np.sum(Weights), ' Gpc-3 yr-1')

        ########################
        # Get the Hist    
        hist, bin_edge = np.histogram(x_vals, weights = Weights, bins=bins)
        y_vals = hist/mass_binw
        center_bins = (bin_edge[:-1] + bin_edge[1:])/2.

        # And the KDE
        kernel = stats.gaussian_kde(x_vals, bw_method=kde_width, weights=Weights)
        binwidth = np.diff(bin_edge)

        m10.append(kernel(10)*sum(hist)) # append value at reference mass 
        m25.append(kernel(25)*sum(hist)) # append value at reference mass 
        m40.append(kernel(40)*sum(hist)) # append value at reference mass 

        ########################
        # Plot the Hist 
        if show_hist:
            plot_lines.append(ax.step(center_bins, y_vals,  where='mid',label = None,#labels[i]+'$%s \mathrm{ \ Gpc^{-3} yr^{-1}}$'%(np.round(np.sum(Weights),1)) , 
                                      alpha=1.0, lw = 3.5, zorder = i, color= Color, 
                                      marker = 'o', markersize = 15) ) #edgecolor=color[i],
            # to prevent overflowing
            min_xkde = min(center_bins[y_vals>5e-4])

        ########################
        # Add KDE
        if show_KDE:
            x_KDE = np.arange(0.1,50.,0.1)
            KDEy_vals =  kernel(x_KDE)*sum(hist) #re-normalize the KDE
            leg_labels.append(labels[nplot]+'$%s$'%(np.round(np.sum(Weights),1)))

            if not show_hist:
                plot_lines.append(ax.plot(x_KDE, KDEy_vals, label = '', color=Color, lw= 5,  zorder =i+1,ls = linestyles[nplot]))

        #     ########################
        #     # Bootstrap   
        #     if bootstrap:
        #         indices = np.arange(len(x_vals))
        #         hist_vals = np.zeros((bootstraps, len(x_KDE)))  #center_bins
        #         for b in progressbar( range(len(hist_vals)), "Bootstrapping "+ labels[i] + ":"):
        #             boot_index = np.random.choice(indices, size=len(indices), replace=True)
        #             kernel         = stats.gaussian_kde(x_vals[boot_index], bw_method=kde_width, weights=Weights[boot_index])
        #             Hist, _        = np.histogram(x_vals[boot_index], bins=bins, weights=Weights[boot_index],density=False)
        #             hist_vals[b]   = kernel(x_KDE)*sum(Hist)

        #         # calculate 1- and 2- sigma percentiles
        #         y_vals = hist_vals#/mass_binw

        #         percentiles = np.percentile(y_vals, [10., 90.], axis=0)
        #         print('nplot',nplot, 'np.shape(percentiles)', np.shape(percentiles))
        #         ax.fill_between(x_KDE, percentiles[0],percentiles[1], alpha=0.4, color=Color, zorder = 11) # 1-sigma

        nplot += 1


    #########################################
    # Show the variation in SFR at 3 different masses
    reference_masses = [10, 25, 40]
    for m, mpoint in enumerate([m10, m25, m40]):
        print('m', np.median(mpoint), max(mpoint), min(mpoint))
        print()
        ax.vlines(x=reference_masses[m], ymin=min(mpoint), ymax=max(mpoint), colors='k', lw=3, zorder = 20)
        ax.hlines(y=[min(mpoint), max(mpoint)], xmin=reference_masses[m]-0.5, xmax=reference_masses[m]+0.5, linewidth=3, color='k', zorder = 20)
        ax.text(reference_masses[m] - 0.7, (max(mpoint)+min(mpoint))/2 , r'%s $\times $'%(np.round( (max(mpoint)/min(mpoint))[0] , 1)) ,
            bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.1', alpha = 0.5), ha = 'right', size = 25, zorder = 20)

    #########################################
    # plot values
    plt.text(0.63, 0.85, titletext, ha = 'center', transform=ax.transAxes, size = 25)

    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)

    #####
    # add legend for simulations
    leg = ax.legend(**leg_args)
    leg = plt.legend([l[0] for l in plot_lines], [l for l in leg_labels ],  **leg_args)
    leg.set_zorder(102)
    leg._legend_box.align = "right"

    # Legend for GWTC-3
    if plot_LIGO:
        plt.gca().add_artist(legend1)

    s = ['$[%s \leq z < %s]$'%(z_bin_edges[a],z_bin_edges[a+1]) for a in range(0,len(z_bin_edges)-1)]   

    ax.set_xlabel(xlabel, fontsize = 30)
    ax.set_ylabel(ylabel, fontsize = 30)

    ax.set_yscale('log')

    if not multipanel:
        if save_plot:
            plt.savefig(save_loc+'/'+save_name , bbox_inches='tight')

        plt.show()
        # clear memory
        gc.collect()
    else:
        return ax



#################################################################################################
#                                                                                               #
#                                          Call plots                                           #
#                                                                                               #
#################################################################################################

fig = plt.figure( figsize = (24, 28))

####################################################
# width of SFRD at z=0
####################################################
#add first subplot in layout that has 3 rows and 2 columns
subplot1 = fig.add_subplot(321)

ax1 = plot_mass_distribution(sim_dir = data_dir, x_key = 'M_moreMassive',  rate_keys  = ['Rates_mu00.025_muz-0.049_alpha-1.778_sigma0%s_sigmaz0.048_a0.017_b1.481_c4.452_d5.913'%(x) for x in [0.7, 1.129, 2.0]],
                       show_hist = False, show_KDE = True, kde_width = 0.1, plot_LIGO = True, Color =  'navy',
                       only_CE = only_CE, only_stable = only_stable, 
                       bootstrap = False, bootstraps = 50, save_name = 'SFRD_width_variations.pdf', titletext = "Width of metallicity dist."+"\n"+r"$\omega_0$, (scale $z=0$)",
                       labels = [r'$\mathrm{Narrow: \ }  (\omega_0 = 0.800) \  \mathcal{R}_{0.2} = \ $',
                                 r'$\mathrm{Fiducial: \ } (\omega_0 = 1.125) \ \mathcal{R}_{0.2}= \ $', 
                                 r'$\mathrm{Wide: \ } \phantom{xx} (\omega_0 = 1.400) \  \mathcal{R}_{0.2} = \ $'],
                      multipanel = True, subplot = subplot1)

# 'Rates_mu00.025_muz-0.049_alpha-3.5_sigma01.129_sigmaz0.048_a0.017_b1.481_c4.452_d5.913
# Rates_mu00.025_muz-0.049_alpha-1.778_sigma01.129_sigmaz0.048_a0.017_b1.481_c4.452_d5.913

####################################################
# Redshift evolution of the width
####################################################
#add Second subplot in layout that has 3 rows and 2 columns
subplot2 = fig.add_subplot(322)

ax2 = plot_mass_distribution(sim_dir = data_dir, x_key = 'M_moreMassive',  rate_keys = ['Rates_mu00.025_muz-0.049_alpha-1.778_sigma01.129_sigmaz%s_a0.017_b1.481_c4.452_d5.913'%(x) for x in [0.0, 0.048, 0.1]],
                       show_hist = False, show_KDE = True, kde_width = 0.1, plot_LIGO = True, Color = '#00a6a0', 
                       only_CE = only_CE, only_stable = only_stable,
                       bootstrap = False, bootstraps = 50, save_name = 'SFRD_zevol_width_variations.pdf',  titletext = "Redshift evol. width of metallicity dist." +"\n"+ r"$\omega_z$, (scale z evol.)",
                       labels = [r'$\mathrm{Flat \ width: \ } \phantom{i} (\omega_z = 0.025) \ \mathcal{R}_{0.2} = \ $',
                                 r'$\mathrm{Fiducial: \ } \phantom{xxi} (\omega_z = 0.050) \ \mathcal{R}_{0.2}= \ $', 
                                 r'$\mathrm{Steep \ width: \ } (\omega_z = 0.100) \ \mathcal{R}_{0.2} = \ $'],
                        multipanel = True, subplot = subplot2)


####################################################
# Mean metallicity at z=0
####################################################
#add third subplot in layout that has 3 rows and 2 columns
subplot3 = fig.add_subplot(323)

ax3 = plot_mass_distribution(sim_dir = data_dir, x_key = 'M_moreMassive',  rate_keys = ['Rates_mu0%s_muz-0.049_alpha-1.778_sigma01.129_sigmaz0.048_a0.017_b1.481_c4.452_d5.913'%(x) for x in [0.007, 0.025, 0.035]],
                       show_hist = False, show_KDE = True, kde_width = 0.1, plot_LIGO = True, Color = '#e1131d', 
                       only_CE = only_CE, only_stable = only_stable,
                       bootstrap = False, bootstraps = 50, save_name = 'SFRD_meanZ_variations.pdf',  titletext = 'Mean metallicity'+"\n"+r"$\mu_0$",
                       labels = [r'$\mathrm{low \ <Z_0> : \ } \phantom{x} (\mu_0 = 0.015) \ \mathcal{R}_{0.2} = \ $',
                                 r'$\mathrm{Fiducial : \ } \phantom{xxxi} (\mu_0 = 0.025) \ \mathcal{R}_{0.2} = \ $',
                                 r'$\mathrm{high \ <Z_0> : \ } \phantom{i} (\mu_0 = 0.035) \ \mathcal{R}_{0.2} = \ $'],
                        multipanel = True, subplot = subplot3)

####################################################
# Redshift evolution of mean metallicity
####################################################
#add 4th subplot in layout that has 3 rows and 2 columns
subplot4 = fig.add_subplot(324)

ax4 = plot_mass_distribution(sim_dir = data_dir, x_key = 'M_moreMassive',  rate_keys = ['Rates_mu00.025_muz%s_alpha-1.778_sigma01.129_sigmaz0.048_a0.017_b1.481_c4.452_d5.913'%(x) for x in [0.0, -0.049, -0.5]],
                       show_hist = False, show_KDE = True, kde_width = 0.1, plot_LIGO = True, Color = '#ff717b', 
                       only_CE = only_CE, only_stable = only_stable,
                       bootstrap = False, bootstraps = 50, save_name = 'SFRD_zevol_mean_variations.pdf', titletext = "Redshift evol. of mean metallicity" +"\n"+ r"$\mu_z$", 
                       labels = [r'$\mathrm{Flat: \ } \phantom{xxi} (\mu_z = -0.01) \ \mathcal{R}_{0.2} = \ $',
                                 r'$\mathrm{Fiducial: \ } (\mu_z = -0.05) \ \mathcal{R}_{0.2}= \ $', 
                                 r'$\mathrm{Steep: \ } \phantom{xx} (\mu_z = -0.25) \ \mathcal{R}_{0.2} = \ $'],
                        multipanel = True, subplot = subplot4)


####################################################
# Skewness
####################################################
#add 5th subplot in layout that has 3 rows and 2 columns
subplot5 = fig.add_subplot(325)

ax5 = plot_mass_distribution(sim_dir = data_dir, x_key = 'M_moreMassive',  rate_keys = ['Rates_mu00.025_muz-0.049_alpha%s_sigma01.129_sigmaz0.048_a0.017_b1.481_c4.452_d5.913'%(x) for x in [0.0, -1.778, -6.0]],
                       show_hist = False, show_KDE = True, kde_width = 0.1, plot_LIGO = True, Color = '#acbf00', 
                       only_CE = only_CE, only_stable = only_stable,
                       bootstrap = False, bootstraps = 50, save_name = 'SFRD_skewness_variations.pdf', titletext = "Skewness of metallicity dist." +"\n"+ r"$\alpha$, (shape)", 
                       labels = [r'$\mathrm{Symmetric: \ } (\alpha = -0.9)   \ \mathcal{R}_{0.2} = \ $',
                                 r'$\mathrm{Fiducial: \  } \phantom{xx} (\alpha = -1.77)  \ \mathcal{R}_{0.2}= \ $', 
                                 r'$\mathrm{Skewed: \    } \phantom{xxi} (\alpha = -3.5)  \ \mathcal{R}_{0.2} = \ $'],
                        multipanel = True, subplot = subplot5)


####################################################
# Star formation norm
####################################################
#add 6th subplot in layout that has 3 rows and 2 columns
subplot6 = fig.add_subplot(326)


ax6 = plot_mass_distribution(sim_dir = data_dir, x_key = 'M_moreMassive',  
                       rate_keys = ['Rates_mu00.025_muz-0.049_alpha-1.778_sigma01.129_sigmaz0.048_a0.01_b2.6_c3.2_d6.2',
                                   'Rates_mu00.025_muz-0.049_alpha-1.778_sigma01.129_sigmaz0.048_a0.017_b1.481_c4.452_d5.913', 
                                   'Rates_mu00.025_muz-0.049_alpha-1.778_sigma01.129_sigmaz0.048_a0.04_b2.5_c2.9_d4.5'],
                       show_hist = False, show_KDE = True, kde_width = 0.1, plot_LIGO = True, Color = '#ecb05b', 
                       only_CE = only_CE, only_stable = only_stable,
                       bootstrap = False, bootstraps = 50, save_name = 'SFRD_skewness_variations.pdf', titletext = "Overall SFR history"+"\n"+ r'$ \mathrm{SFRD(}z\rm{)} \ [a,b,c,d]$', 
                       labels = [r'$\mathrm{Madau \ \& \ Fragos \ 2017: } \phantom{xxx} \ \mathcal{R}_{0.2}= \ $', 
                                 r'$\mathrm{Fiducial: \ } \phantom{xxxxxxxxx} \ \mathcal{R}_{0.2}= \ $', 
                                 r'$\mathrm{Max \ SB: \ B18/C17, \ Chruslinska \ et \ al. \ 2021:}  \ \mathcal{R}_{0.2} = \ $'],
                        multipanel = True, subplot = subplot6)



####################################################
# Final plot properties
fig.savefig(save_loc + '/Mass_distributions_'+ channel_string+'_SFRD_variations.pdf' , bbox_inches='tight')
# plt.show()


