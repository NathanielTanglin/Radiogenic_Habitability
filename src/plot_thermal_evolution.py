from lib import *
from paths import path

#planet_dynamics = read_vplanet(path("src", "parameter_sweeps", "kv_atm_central_hz", "kv.earth.forward"))
planet_dynamics = read_vplanet(path("data", "parameter_sweeps", "synthesis", "kv_atm", "run_K1_ThU1_a1", "star.earth.forward"))
save_path = path("plots", "kv_atm_outer_hz_hot.png")
plot_thermal_evolution(save_path, planet_dynamics = planet_dynamics)

#planet_dynamics = read_vplanet(os.path.join('parameter_sweeps', 'sun_atm_1_au_slow', 'sol.earth.forward'))
#save_path = os.path.join(os.pardir, 'plots', 'sun_atm_1_au_slow.png')
#plot_thermal_evolution(save_path, planet_dynamics = planet_dynamics)