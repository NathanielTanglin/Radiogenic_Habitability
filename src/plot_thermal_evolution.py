from lib import *

planet_dynamics = read_vplanet(os.path.join('parameter_sweeps', 'sun_atm_1_au_highO', 'sol.earth.forward'))
save_path = os.path.join(os.pardir, 'plots', 'sun_atm_1_au_highO.png')
plot_thermal_evolution(save_path, planet_dynamics = planet_dynamics)

planet_dynamics = read_vplanet(os.path.join('parameter_sweeps', 'sun_atm_1_au_slow', 'sol.earth.forward'))
save_path = os.path.join(os.pardir, 'plots', 'sun_atm_1_au_slow.png')
plot_thermal_evolution(save_path, planet_dynamics = planet_dynamics)