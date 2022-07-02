"""
This script makes a slrum job to run FastCosmicIntegration with a specified set of 
"""
import numpy as np
import os
from subprocess import Popen, PIPE, call
import subprocess
import sys
import paths
import time


#################################################################
## 
##   Should be Changed by user ##
##
#################################################################
#################################################################
def init():
    # define initial settings
    global data_dir, script_dir, COMPASfilename, rate_file_name, user_email, fid_dpdZ_parameters, fid_sfr_parameters, SlurmJobString
    global mu0_best, muz_best, sigma0_best, sigmaz_best, alpha0_best,sf_a_best, sf_b_best, sf_c_best, sf_d_best
    # there is an extra /src in the paths. due to the import location (remove it)
    data_dir   = str(paths.data)[0:-8] + 'data/'
    script_dir = str(paths.scripts)[0:-11] + 'scripts/'

    COMPASfilename  = 'small_COMPAS_Output_wWeights.h5'#'COMPAS_Output_wWeights.h5'
    rate_file_name  = 'small_Rate_info.hdf5'#'Rate_info.hdf5'
    user_email      = "aac.van.son@gmail.com"


    mu0_best, muz_best, sigma0_best, sigmaz_best, alpha0_best,sf_a_best, sf_b_best, sf_c_best, sf_d_best = np.loadtxt(data_dir+ '/best_fit_parameters.txt',unpack=True)

    fid_dpdZ_parameters = [mu0_best, muz_best, sigma0_best, sigmaz_best, alpha0_best]
    fid_sfr_parameters  = [sf_a_best, sf_b_best, sf_c_best, sf_d_best]


    ##################################################################
    # This is the slurm script youre using
    #SBATCH --partition=%s              # Partition to submit to
    ##################################################################
    # note indentation needs to be like this
    SlurmJobString="""#!/bin/bash
#SBATCH --job-name=%s          #job name
#SBATCH --nodes=%s             # Number of nodes
#SBATCH --ntasks=%s            # Number of cores
#SBATCH --output=%s            # output storage file
#SBATCH --error=%s             # error storage file
#SBATCH --time=%s              # Runtime in minutes
#SBATCH --mem=%s               # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH -p %s
#SBATCH --mail-user=%s         # Send email to user
#SBATCH --mail-type=FAIL       #
#
#Print some stuff on screen
echo $SLURM_JOB_ID
echo $SLURM_JOB_NAME
echo $SLURM_ARRAY_TASK_ID
#
#Set variables
export QT_QPA_PLATFORM=offscreen # To avoid the X Display error
#
cd %s
#
# Run your job
%s
    """

###############################################
###
###############################################
def RunSlurmBatch(run_dir = None, job_name = "stroopwafel_interface", dependency = False, dependent_ID = None):
    if not dependency:
        sbatchArrayCommand = 'sbatch ' + os.path.join(run_dir+job_name+'.sbatch') 
    else:
        sbatchArrayCommand = 'sbatch --dependency=afterok:' + str(int(dependent_ID)) + ' ' + os.path.join(run_dir+job_name+'.sbatch') 

    # Open a pipe to the sbatch command.
    proc = Popen(sbatchArrayCommand, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)

    # Send job_string to sbatch
    if (sys.version_info > (3, 0)):
        proc.stdin.write(sbatchArrayCommand.encode('utf-8'))
    else:
        proc.stdin.write(sbatchArrayCommand)

    print('sbatchArrayCommand:', sbatchArrayCommand)
    out, err = proc.communicate()
    print("out = ", out)
    job_id = out.split()[-1]
    print("job_id", job_id)
    return job_id


#########################################################
# Assuming you have The COMPAS simulation data output
# Make a slurm job to call 
def Call_Cosmic_Integration(root_out_dir, COMPASfilename, rate_file_name, jname = 'CI',
                            ZdepSFRD_param_sets = [[],[]], #[[fid_dpdZ_parameters, fid_sfr_parameters]],
                           partitions = 'demink,conroy,hernquist,shared', Wtime = "5:00:00", mem = "150000",
                           number_of_nodes = 1, number_of_cores = 1):
    """
    Call slurm batch 

    Args:
        root_out_dir           --> [string] Path to the COMPAS file that contains the simulation data
        COMPASfilename         --> [string] Name of the COMPAS file
        rate_file_name         --> [string] Name of the output file containing the rates
        ZdepSFRD_param_sets    --> [array of float arrays] 
        # consists of: 
            dco_type               --> [string] Which DCO type you used to calculate rates 
            mu0                    --> [float]  metallicity dist: expected value at redshift 0
            muz                    --> [float]  metallicity dist: redshift evolution of expected value
            sigma0                 --> [float]  metallicity dist: width at redshhift 0
            sigmaz                 --> [float]  metallicity dist: redshift evolution of width
            alpha                  --> [float]  metallicity dist: skewness (0 = lognormal)

        maxzdet                   --> [float] max z for detection rate, will be used as maxz: Maximum redshhift up to where we would like to store the data

    """
    # Index to name your job
    n_CI = 1 

    # Run over each complete variation of the metallicity dependent SFRD
    for SFRD_zZ in ZdepSFRD_param_sets:  
        job_name = jname+str(n_CI)
        
        print(10* "*" + ' You are Going to Run FastCosmicIntegration.py')
        print('dPdZ data = ', SFRD_zZ[0])
        print('SFR(z) data = ', SFRD_zZ[1])
        
        mu0, muz, sigma0, sigmaz, alpha0 = SFRD_zZ[0]
        sf_a, sf_b, sf_c, sf_d           = SFRD_zZ[1]
        DEPEND, append_job_id = False, 0

        # Flag to pass to FasCosmicIntegrator
        Flags = " --path "+root_out_dir + " --filename "+COMPASfilename+" --outfname " +rate_file_name+\
        " --mu0 " +str(mu0)+" --muz "+str(muz)+" --sigma0 "+str(sigma0)+" --sigmaz "+str(sigmaz)+" --alpha "+str(alpha0)+\
        " --aSF " +str(sf_a)+" --bSF "+str(sf_b)+" --cSF "+str(sf_c)+" --dSF "+str(sf_d)+\
        " --weight "+"mixture_weight"+ " --zstep "+"0.01"+" --sens "+"O3"+ " --m1min "+"10."+ " --dco_type BBH"+\
        " --redshiftBinSize "+"0.05" + ' --maxzdet ' + "0.5" #+ " --BinAppend "

        run_dir = script_dir +'/CosmicIntegration/'

        # Make and safe a slurm command
        job_line = "python FastCosmicIntegration.py "+Flags+" > "+ root_out_dir + "/slurm_out/"+job_name+".log"

        # Make slurm script string
        interface_job_string = SlurmJobString % (job_name, number_of_nodes, number_of_cores, \
        root_out_dir+'/slurm_out/'+job_name+'.out', root_out_dir+'/slurm_out/'+job_name+'.err', Wtime, mem, partitions, user_email, run_dir, job_line)
        
        # Write your bash file
        sbatchFile = open(run_dir+job_name+'.sbatch','w')
        print('writing ', run_dir+job_name+'.sbatch')
        sbatchFile.write(interface_job_string)
        sbatchFile.close()
  
        # Submit the job to sbatch! 
        CIjob_id = RunSlurmBatch(run_dir = run_dir, job_name = job_name ,\
        dependency = DEPEND, dependent_ID = append_job_id)

        n_CI += 1
        DEPEND, append_job_id = False, CIjob_id # no need for it to be dependent
            

    ###########################
    # Now wait for your (last) job to be done
    command = "sacct  -j %s --format State "%(CIjob_id.decode("utf-8"))
    print(command)
    done = False
    while not done:
        # First check the status of your job with sacct
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        nline = 0
        while True:
            line = p.stdout.readline()
            nline +=1
            #print(line)
            if not line:
                break
            if nline == 3:
                break


        result = str(line)
        print('result = ', result)
        if b"COMPLETE" in line:
            print('YAY your job finished!')
            done = True
        elif b"FAILED" in line:
            print('Job failed :( %s')%CIjob_id
            done = True
        elif np.logical_or(b"RUNNING" in line, b"PENDING" in line):
            print('darn, still running, check back in 10 sec')
            time.sleep(150) # Sleep 2.5 min and then check back


    print(10* "*" + " You are all done with this job! " + 10* "*")
    
    


#################################################################
# All at once
#################################################################
#Call_Cosmic_Integration(data_dir, COMPASfilename, rate_file_name, 
#                        ZdepSFRD_param_sets =[[fid_dpdZ_parameters, fid_sfr_parameters],
#                                              [[0.015, muz_best, sigma0_best, sigmaz_best, alpha0_best], fid_sfr_parameters], # mu0_variations
#                                             [[0.035, muz_best, sigma0_best, sigmaz_best, alpha0_best], fid_sfr_parameters],
#                                              [[mu0_best, -0.01, sigma0_best, sigmaz_best, alpha0_best], fid_sfr_parameters],# muz_variations
#                                             [[mu0_best, -0.25, sigma0_best, sigmaz_best, alpha0_best], fid_sfr_parameters],
#                                              [[mu0_best, muz_best, 0.8, sigmaz_best, alpha0_best], fid_sfr_parameters],# omega0_variations
#                                             [[mu0_best, muz_best, 1.4, sigmaz_best, alpha0_best], fid_sfr_parameters],
#                                              [[mu0_best, muz_best, sigma0_best, 0.025, alpha0_best], fid_sfr_parameters],# omegaz_variations
#                                             [[mu0_best, muz_best, sigma0_best, 0.1, alpha0_best], fid_sfr_parameters],
#                                              [[mu0_best, muz_best, sigma0_best, sigmaz_best, -0.9], fid_sfr_parameters],# alpha_variations
#                                             [[mu0_best, muz_best, sigma0_best, sigmaz_best, -3.5], fid_sfr_parameters],
#                                              [fid_dpdZ_parameters, [0.01, 2.60, 3.20, 6.20]],# SFR(z) variations
#                                              [fid_dpdZ_parameters, [0.01, 2.77, 2.90, 4.70]] ],
#                        partitions = 'demink,conroy,hernquist,shared', Wtime = "1:00:00", mem = "120000")

