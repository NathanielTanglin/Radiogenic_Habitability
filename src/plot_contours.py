from lib import *
from paths import *

#build_contour_VPLanet(save_name = path('plots', 'earth_sun_inner_hz.png'), directory = path('data', 'parameter_sweeps', 'Earth_Sun_Inner_HZ'), planet_file_name = 'sol.earth.forward')

# Same errors. Indexing problem. Related to choice of binary flags (e.g. DESSICATION_TIME).
#build_contour_VPLanet(save_name = path('plots', 'earth_sun_outer_hz.png'), directory = path('data', 'parameter_sweeps', 'Earth_Sun_Outer_HZ'), planet_file_name = 'sol.earth.forward')

#build_contour_VPLanet(save_name = path('plots', 'earth_trappist_central_hz.png'), directory = path('data', 'parameter_sweeps', 'Earth_Trappist_Central_HZ'), planet_file_name = 'trappist.earth.forward', mode = ATM | DESSIC_TIME)
#build_contour_VPLanet(save_name = path('plots', 'earth_trappist_central_hz_water.png'), directory = path('data', 'parameter_sweeps', 'Earth_Trappist_Central_HZ_Water'), planet_file_name = 'trappist.earth.forward', mode = WATER | DESSIC_TIME)

#build_contour_VPLanet(save_name = path('plots', 'earth_sun_1_AU.png'), directory = path('data', 'parameter_sweeps', 'Earth_Sun_1_AU'), planet_file_name = 'sol.earth.forward', mode = ATM | PERCENT)
#build_contour_VPLanet(save_name = path('plots', 'earth_sun_1_AU_water.png'), directory = path('data', 'parameter_sweeps', 'Earth_Sun_1_AU_Water'), planet_file_name = 'sol.earth.forward', mode = WATER)

#build_contour_VPLanet(save_name = path('plots', 'kv_atm_central_hz.png'), directory = path('data', 'parameter_sweeps', 'kv_atm_central_hz'), planet_file_name = 'kv.earth.forward', mode = ATM | DESSIC_TIME)
build_contour_VPLanet(save_name = path('plots', 'kv_water_central_hz.png'), directory = path('data', 'parameter_sweeps', 'kv_water_central_hz'), planet_file_name = 'kv.earth.forward', mode = WATER)