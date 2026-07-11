import os
from paths import path

GKM_host_parameter_sweep = path("src", "run_GKM_host_parameter_sweeps.ipynb")
os.system("python " + GKM_host_parameter_sweep)

# Runs the python script for the parameter sweep which varies stellar host mass.
stellar_host_parameter_sweep = path("src", "run_stellar_host_parameter_sweeps.py")
os.system("python " + stellar_host_parameter_sweep)