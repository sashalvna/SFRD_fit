"""
This script makes a slrum job to run FastCosmicIntegration with a specified set of 
"""
#####################################
# dpdZ_parameters = [mu0_best, muz_best, sigma0_best, sigmaz_best, alpha0_best]
#####################################
mu0_best     = 0.025
muz_best     = -0.049
sigma0_best  = 1.129
sigmaz_best  = 0.048
alpha0_best  = -1.778

#####################################
# sfr_parameters  = [sf_a_best, sf_b_best, sf_c_best, sf_d_best]
#####################################
sf_a_best     = 0.017
sf_b_best     = 1.481
sf_c_best     = 4.452
sf_d_best     = 5.913

user_email = 'aac.van.son@gmail.com'

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
def Call_Cosmic_Integration(root_out_dir, COMPASfilename, rate_file_name, 
                            ZdepSFRD_param_sets = [[dpdZ_parameters, sfr_parameters]]
                           partitions = 'demink,conroy,hernquist,shared', Wtime = "5:00:00", mem = "150000"):
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
                
            maxz                   --> [float] Maximum redshhift up to where we would like to store the data
            sensitivity            --> [string] Which detector sensitivity you used to calculate rates 
            append_binned_by_z     --> [Bool] to save space, bin rates by redshiftbin and append binned rates
            redshift_binsize       --> [float] if append_binned_by_z, how big should your redshift bin be

        Returns:
            h_new                  --> [hdf5 file] Compas output file with a new group "rates" with the same shape as DoubleCompactObjects x redshifts
    """
    
    for SFRD_zZ in ZdepSFRD_param_sets:
        mu0, muz, sigma0, sigmaz, alpha0 = SFRD_zZ[0]
        sf_a, sf_b, sf_c, sf_d           = SFRD_zZ[1]
    
        print(10* "*" + ' You are Going to Run FastCosmicIntegration.py')
        DEPEND, append_job_id = False, 0

        # Run over your metallicity density parameters of interest
        n_CI = 1

        Flags = " --path "+root_out_dir+"/output/"+" --filename "+COMPASfilename+\
        " --mu0 " +str(mu0)+" --muz "+str(muz)+" --sigma0 "+str(sig0)+" --sigmaz "+str(sigz)+" --alpha "+str(al)+\
        " --aSF " +str(sf_a)+" --bSF "+str(sf_b)+" --cSF "+str(sf_c)+" --dSF "+str(sf_d)+\
        " --weight "+"mixture_weight"+ " --zstep "+"0.01"+" --sens "+"O3"+ " --m1min "+"10."+ " --dco_type BBH"+\
        " --BinAppend "+ " --redshiftBinSize "+"0.05" + '--maxzdet' + "0.5"

        if cluster_run:
            # Make and safe a slurm command
            CI_job_string = MakeSlurmBatch(OUT_DIR = root_out_dir, sub_dir = 'CosmicIntegration/', python_name = "FastCosmicIntegration",\
            job_name = "COMPAS_CI"+str(n_CI), number_of_nodes = 1, number_of_cores = 1, partition=partitions,\
            walltime = Wtime, memory = mem, email = user_email, flags= Flags)

            # Submit the job to sbatch! 
            CIjob_id = RunSlurmBatch(run_dir = root_out_dir+'/masterfolder/CosmicIntegration/', job_name = "COMPAS_CI"+str(n_CI),\
            dependency = DEPEND, dependent_ID = append_job_id)
            n_CI += 1
            DEPEND, append_job_id = True, CIjob_id
            
        else:
            print('You will run the CI locally')
            # Change the current working directory
            os.chdir(root_out_dir+'masterfolder/CosmicIntegration')
            # execute this job locally (slow)
            job_line = "python FastCosmicIntegration.py "+Flags+" > "+"COMPAS_CI"+str(n_CI)+".log"
            print('job_line', job_line)

            with open("./COMPAS_CI"+str(n_CI)+".err", "w+") as f:
                subprocess.call(job_line, shell=True, stdin=PIPE, stdout=f, stderr=f)


    print(10* "*" + " You are all done with this job! " + 10* "*")
    
##################################################################
    
    
    
    
    
    