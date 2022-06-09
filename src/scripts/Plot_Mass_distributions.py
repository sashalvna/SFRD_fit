"""

# Plotting the BH mass distribution for several SFRD Z-distribution variations

"""

import sys
import numpy as np
import matplotlib.pyplot as plt

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


######################################
## locations
save_loc    =  '/Users/lieke/surfdrive/Documents'+'/SFRD_fit/src/tex/figures/' #/n/home04/lvanson/
data_dir    = '/Volumes/StorageSpac'+'/CompasOutput/v02.19.04/SFRD_fit_data/fWR1.0coolWind1.0/output/' # '/n/holystore01/LABS/hernquist_lab/Users/lvanson/'


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
                   bootstrap = False, bootstraps = 10, x_lim=(0.,50),  y_lim = (1e-2,50), 
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
		LIGO_data_dir = '/Volumes/StorageSpac/CompasOutput//output/'#'/n/home04/lvanson/LowMBH_peak/output/'
		#################################################
		## prep to grab Powerlaw + Peak data from O3
		#################################################
		input_fname = LIGO_data_dir+'/GWTC-3-population-data/analyses/PowerLawPeak/o3only_mass_c_iid_mag_iid_tilt_powerlaw_redshift_mass_data.h5'
		output_fname = './plots'
		cli_parser = argparse.ArgumentParser()
		cli_parser.add_argument("input_fname")
		cli_parser.add_argument("output_fname")
		cli_args = cli_parser.parse_args( (input_fname + ' ' + output_fname).split() )
		O3_only_PPD = cli_args.input_fname
		O3_only_result = O3_only_PPD.replace("_mass_data.h5", "_result.json")

		with open(O3_only_result, "r") as jfile:
			plpeak_jf = json.load(jfile)

		##############################
		# Calculate relevant values
		mass_1 = np.linspace(2, 100, 1000)
		mass_ratio = np.linspace(0.1, 1, 500)
		color_plpeak = 'grey'#'#1f78b4'

		with h5.File(O3_only_PPD, "r") as f:
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
	DCO = mfunc.read_data(loc = sim_dir +'/COMPAS_Output_wWeights.h5')

	# We'll show the change in rate between SFRD at several reference masses
	m10, m25, m40 = [], [], []

	print('nplot', nplot, '\n')
	####################################################
	### Loop over SFRD
	for i, rate_key in enumerate(rate_keys):
		print('rate_key', rate_key)

		# ### ## Reading Rate data ##
		# try:
		#     DCO_mask, redshifts, intrinsic_rate_density, intrinsic_rate_density_z0  = mfunc.read_rate_data(loc = sim_dir + '/Rate_info.hdf5', rate_key = rate_key)
		## Open hdf5 file
		with h5.File(sim_dir + '/Rate_info.hdf5' ,'r') as File:
			redshifts                 = File['redshifts'][()]
			# Different per rate key:
			DCO_mask                  = File[rate_key]['DCOmask'][()] # Mask from DCO to merging systems  
			#(contains filter for RLOF>CE and optimistic CE)
			intrinsic_rate_density    = File[rate_key]['merger_rate'][()]

		# except:
			# print('\n error reading', rate_key)
			# continue
	    
		# # # # # # # # # # # # # # # # # # 
		#first bring it to the same shape as the rate table
		merging_BBH    = DCO[DCO_mask]
		#then apply the additional mask based on your prefs
		merging_BBH         = merging_BBH[(merging_BBH['Stellar_Type@ZAMS(1)'] != 16)]
		Red_intr_rate_dens  = intrinsic_rate_density[(DCO['Stellar_Type@ZAMS(1)'][DCO_mask] != 16), :]

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
	plt.text(0.90, 0.90, titletext, ha = 'right', transform=ax.transAxes, size = 25)

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
#    																							#
#       Call plots																				#	
#																								#
#################################################################################################

fig = plt.figure( figsize = (24, 30))


####################################################
# width of SFRD at z=0
####################################################
#add first subplot in layout that has 3 rows and 2 columns
subplot1 = fig.add_subplot(321)

ax1 = plot_mass_distribution(sim_dir = data_dir, x_key = 'M_moreMassive',  rate_keys = ['Rates_mu00.025_muz-0.05_alpha-1.77_sigma0%s_sigmaz0.05_zBinned'%(x) for x in [0.8, 1.125, 1.4]],
                       show_hist = False, show_KDE = True, kde_width = 0.07, plot_LIGO = True, Color =  'navy', 
                       bootstrap = False, bootstraps = 50, save_name = 'SFRD_width_variations.pdf', titletext = 'Width of metallicity  dist. z=0',
                       labels = [r'$\mathrm{Narrow: \ } \phantom{xxx} (\omega_0 = 0.800) \  \mathcal{R}_{0} = \ $',
                                 r'$\mathrm{Fiducial: \ } \phantom{xxi} (\omega_0 = 1.125) \ \mathcal{R}_{0}= \ $', 
                                 r'$\mathrm{Wide \ SFRD: \ }  (\omega_0 = 1.400) \  \mathcal{R}_{0} = \ $'],
                      multipanel = True, subplot = subplot1)


####################################################
# Redshift evolution of the width
####################################################
#add Second subplot in layout that has 3 rows and 2 columns
subplot2 = fig.add_subplot(322)

ax2 = plot_mass_distribution(sim_dir = data_dir, x_key = 'M_moreMassive',  rate_keys = ['Rates_mu00.025_muz-0.05_alpha-1.77_sigma01.125_sigmaz%s_zBinned'%(x) for x in [0.025, 0.05, 0.1]],
                       show_hist = False, show_KDE = True, kde_width = 0.07, plot_LIGO = True, Color = '#00a6a0', 
                       bootstrap = False, bootstraps = 50, save_name = 'SFRD_zevol_width_variations.pdf',  titletext = 'z-evol  of  metallicity  width',
                       labels = [r'$\mathrm{Flat \ width: \ } \phantom{i} (\omega_z = 0.025) \ \mathcal{R}_{0} = \ $',
                                 r'$\mathrm{Fiducial: \ } \phantom{xxi} (\omega_z = 0.050) \ \mathcal{R}_{0}= \ $', 
                                 r'$\mathrm{Steep \ width: \ } (\omega_z = 0.100) \ \mathcal{R}_{0} = \ $'],
                        multipanel = True, subplot = subplot2)


####################################################
# Mean metallicity at z=0
####################################################
#add third subplot in layout that has 3 rows and 2 columns
subplot3 = fig.add_subplot(323)

ax3 = plot_mass_distribution(sim_dir = data_dir, x_key = 'M_moreMassive',  rate_keys = ['Rates_mu0%s_muz-0.05_alpha-1.77_sigma01.125_sigmaz0.05_zBinned'%(x) for x in [0.015, 0.025, 0.035]],
                       show_hist = False, show_KDE = True, kde_width = 0.07, plot_LIGO = True, Color = '#e1131d', 
                       bootstrap = False, bootstraps = 50, save_name = 'SFRD_meanZ_variations.pdf',  titletext = 'Mean metallicity z=0',
                       labels = [r'$\mathrm{low \ <Z_0> : \ } \phantom{x} (\mu_0 = 0.015) \ \mathcal{R}_{0} = \ $',
                       			 r'$\mathrm{Fiducial : \ } \phantom{xxxi} (\mu_0 = 0.025) \ \mathcal{R}_{0} = \ $',
                                 r'$\mathrm{high \ <Z_0> : \ } \phantom{i} (\mu_0 = 0.035) \ \mathcal{R}_{0} = \ $'],
                        multipanel = True, subplot = subplot3)


####################################################
# Redshift evolution of mean metallicity
####################################################
#add 4th subplot in layout that has 3 rows and 2 columns
subplot4 = fig.add_subplot(324)

ax4 = plot_mass_distribution(sim_dir = data_dir, x_key = 'M_moreMassive',  rate_keys = ['Rates_mu00.025_muz%s_alpha-1.77_sigma01.125_sigmaz0.05_zBinned'%(x) for x in [-0.01, -0.05, -0.25]],
                       show_hist = False, show_KDE = True, kde_width = 0.07, plot_LIGO = True, Color = '#ff717b', 
                       bootstrap = False, bootstraps = 50, save_name = 'SFRD_zevol_mean_variations.pdf', titletext = 'z-evol of mean metallicity', 
                       labels = [r'$\mathrm{Flat: \ } \phantom{xxi} (\mu_z = -0.01) \ \mathcal{R}_{0} = \ $',
                                 r'$\mathrm{Fiducial: \ } (\mu_z = -0.05) \ \mathcal{R}_{0}= \ $', 
                                 r'$\mathrm{Steep: \ } \phantom{xx} (\mu_z = -0.25) \ \mathcal{R}_{0} = \ $'],
                        multipanel = True, subplot = subplot4)


####################################################
# Skewness
####################################################
#add 4th subplot in layout that has 3 rows and 2 columns
subplot5 = fig.add_subplot(325)

ax5 = plot_mass_distribution(sim_dir = data_dir, x_key = 'M_moreMassive',  rate_keys = ['Rates_mu00.025_muz-0.05_alpha%s_sigma01.125_sigmaz0.05_zBinned'%(x) for x in [-0.9, -1.77, -3.5]],
                       show_hist = False, show_KDE = True, kde_width = 0.07, plot_LIGO = True, Color = '#acbf00', 
                       bootstrap = False, bootstraps = 50, save_name = 'SFRD_skewness_variations.pdf', titletext = 'Skewness of metallicity dist.', 
                       labels = [r'$\mathrm{Symmetric: \ } (\alpha = -0.9)   \ \mathcal{R}_{0} = \ $',
                                 r'$\mathrm{Fiducial: \  } \phantom{xx} (\alpha = -1.77)  \ \mathcal{R}_{0}= \ $', 
                                 r'$\mathrm{Skewed: \    } \phantom{xxi} (\alpha = -3.5)  \ \mathcal{R}_{0} = \ $'],
                        multipanel = True, subplot = subplot5)


####################################################
# Star formation norm
####################################################
#add 4th subplot in layout that has 3 rows and 2 columns
subplot6 = fig.add_subplot(326)


ax6 = plot_mass_distribution(sim_dir = data_dir, x_key = 'M_moreMassive',  
                       rate_keys = ['Rates_mu00.025_muz-0.05_alpha-1.77_sigma01.125_sigmaz0.05_a0.01_b2.6_c3.2_d6.2_zBinned',
                                    'Rates_mu00.025_muz-0.05_alpha-1.77_sigma01.125_sigmaz0.05_zBinned',
                                    'Rates_mu00.025_muz-0.05_alpha-1.77_sigma01.125_sigmaz0.05_a0.01_b2.77_c2.9_d4.7_zBinned'],
                       show_hist = False, show_KDE = True, kde_width = 0.07, plot_LIGO = True, Color = '#ecb05b', 
                       bootstrap = False, bootstraps = 50, save_name = 'SFRD_skewness_variations.pdf', titletext = 'Magnitude of SFR(z)', 
                       labels = [r'$\mathrm{Madau \ \& \ Fragos \ 2017: } \ \mathcal{R}_{0}= \ $', 
                                 r'$\mathrm{Fiducial: \ } \phantom{xxxxxxxxx} \ \mathcal{R}_{0}= \ $', 
                                 r'$\mathrm{Neijssel \ et \ al. \ 2019:  \phantom{xi}  }  \ \mathcal{R}_{0} = \ $'],
                        multipanel = True, subplot = subplot6)


####################################################
# Final plot properties
fig.savefig(save_loc+'/SFRD_variations_combined.pdf' , bbox_inches='tight')
# plt.show()


