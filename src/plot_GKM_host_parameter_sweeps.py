from lib import build_contour_VPLanet
from paths import path

build_contour_VPLanet(save_name = path('plots', 'trappist_atm_central_hz.png'), directory = path('data', 'parameter_sweeps', 'trappist_atm_central_hz'), planet_file_name = 'trappist.earth.forward', mode = ATM | DESSIC_TIME)
build_contour_VPLanet(save_name = path('plots', 'trappist_water_central_hz.png'), directory = path('data', 'parameter_sweeps', 'trappist_water_central_hz'), planet_file_name = 'trappist.earth.forward', mode = WATER | DESSIC_TIME)

build_contour_VPLanet(save_name = path('plots', 'sun_atm_1_au.png'), directory = path('data', 'parameter_sweeps', 'sun_atm_1_au'), planet_file_name = 'sol.earth.forward', mode = ATM)
build_contour_VPLanet(save_name = path('plots', 'sun_water_1_au.png'), directory = path('data', 'parameter_sweeps', 'sun_water_1_au'), planet_file_name = 'sol.earth.forward', mode = WATER)

build_contour_VPLanet(save_name = path('plots', 'kv_atm_central_hz.png'), directory = path('data', 'parameter_sweeps', 'kv_atm_central_hz'), planet_file_name = 'kv.earth.forward', mode = ATM)
build_contour_VPLanet(save_name = path('plots', 'kv_water_central_hz.png'), directory = path('data', 'parameter_sweeps', 'kv_water_central_hz'), planet_file_name = 'kv.earth.forward', mode = WATER)
