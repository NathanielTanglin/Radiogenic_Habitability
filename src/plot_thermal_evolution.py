from lib import *
from paths import path
import os

simulation_directory = path("src", "parameter_sweeps", "sun_atm_1_au")
output_file = "sol.earth.forward"

# Runs the vplanet simulation if necessary.
if output_file not in os.listdir(simulation_directory):
    os.chdir(simulation_directory)
    os.system("vplanet vpl.in")

output_file = path(simulation_directory, output_file)
planet_dynamics = read_vplanet(output_file)

save_path = path("plots", "fiducial_evolution.png")
plot_thermal_evolution(save_path, planet_dynamics = planet_dynamics)