import numpy as np
import os
import copy
from parameter_sweep import Parameter_Sweep
from paths import path

num_cores = ''

while not num_cores.isdigit():
    if num_cores != '':
        print('Please input a valid number!')

    num_cores = input('How many cores would you like to run the parameter sweeps on? Cores: ')

def run_multiplanet():
    os.chdir(path("data", "parameter_sweeps", "synthesis"))
    os.system("pwd")
    os.system('multiplanet -c {cores} -v vspace.in'.format(cores = num_cores))

sweeps = dict()

def add_parameter_sweep(stellar_host_label, hz_range):
    atm_run_directory = stellar_host_label + "_atm"
    water_run_directory = stellar_host_label + "_water"
    
    atm_run_path = path('src', 'parameter_sweeps', 'synthesis', atm_run_directory)
    water_run_path = path('src', 'parameter_sweeps', 'synthesis', water_run_directory)

    sweep_atm = {
        'trial_name': 'run',

        "paths": [path(atm_run_path, 'earth.in'), path(atm_run_path, stellar_host_label + ".in"), path(atm_run_path, 'vpl.in')],

        'earth.in': {
            'input_options': [
                ['-d40KNumMan', '-d40KNumCore'], # Negatives indicate default units for these parameters.
                ['-d232ThNumMan', '-d235UNumMan', '-d238UNumMan'],
                ['dSemi']
            ],
            'names': [
                'K',
                'ThU',
                'a',
            ],
            'ranges': [
                np.linspace(0.0348016632215, 2.6302679919, 2),
                np.linspace(0.68438423992, 5.6767530243, 2),
                hz_range # (inner, outer) AU (fill in on copy)
            ]
        },
    }
    
    sweep_water = copy.deepcopy(sweep_atm)
    sweep_water['paths'] = [path(water_run_path, 'earth.in'), path(water_run_path, stellar_host_label + ".in"), path(water_run_path, 'vpl.in')]

    sweeps[atm_run_directory] = sweep_atm
    sweeps[water_run_directory] = sweep_water

add_parameter_sweep("m1p2", (1.550, 2.733))
add_parameter_sweep("m1p1", (1.287, 2.269))
add_parameter_sweep("sun", (0.9553, 1.685))
add_parameter_sweep("m0p9", (0.6995, 1.234))
add_parameter_sweep("kv", (0.6076, 1.071))
add_parameter_sweep("m0p8", (0.5103, 0.8999))
add_parameter_sweep("m0p7", (0.3659, 0.6453))
add_parameter_sweep("m0p6", (0.2568, 0.4529))
add_parameter_sweep("m0p5", (0.1816, 0.3203))
add_parameter_sweep("m0p4", (0.1341, 0.2365))
add_parameter_sweep("m0p3", (0.09967, 0.1757))
add_parameter_sweep("m0p2", (0.06620, 0.1167))
add_parameter_sweep("m0p1", (0.02829, 0.04989))
add_parameter_sweep("trappist", (0.02237, 0.03945))

def run_sweeps():
    for (sweep_dir, sweep_parameters) in sweeps.items():
        Parameter_Sweep(**sweep_parameters).generate_input_files(path("data", "parameter_sweeps", "synthesis", sweep_dir))
        run_multiplanet()

run_sweeps()