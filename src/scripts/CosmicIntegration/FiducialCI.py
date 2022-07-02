"""
Run a variation of the SFRD 
(this is in a separate file so you dont have to rerun everything if one job fails)
"""
import numpy as np
import os
from subprocess import Popen, PIPE, call
import subprocess
import sys
import paths

import CallCosmicIntegration as CI
from CallCosmicIntegration import init


#################################################################
# fiducial
CI.Call_Cosmic_Integration(init.data_dir, init.COMPASfilename, init.rate_file_name, jname = 'fid',
                        ZdepSFRD_param_sets =[[init.fid_dpdZ_parameters, init.fid_sfr_parameters]],
                        partitions = 'demink,conroy,hernquist,shared', Wtime = "1:00:00", mem = "120000")
                                               
import time
time.sleep(3000) # Sleep until the coscmic integration slurms should be done



