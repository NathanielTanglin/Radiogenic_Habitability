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

add_parameter_sweep("trappist", (0.024967375866069036, 0.050430828529552864))
add_parameter_sweep("m0p1", (0.031060592538580522, 0.06229305781670558))
add_parameter_sweep("m0p2", (0.07235421109392776, 0.1421672474441826))
add_parameter_sweep("m0p3", (0.10877635393888095, 0.2121318467235588))
add_parameter_sweep("m0p4", (0.14612964266104495, 0.28349597780127755))
add_parameter_sweep("m0p5", (0.19761127501413284, 0.3802637544815201))
add_parameter_sweep("m0p6", (0.2776971625930472, 0.5262273650739901))
add_parameter_sweep("m0p7", (0.39124504606264093, 0.7254171525991292))
add_parameter_sweep("m0p8", (0.5359642927048134, 0.9735022355856436))
add_parameter_sweep("kv", (0.6401830775621882, 1.1510924277633037))
add_parameter_sweep("m0p9", (0.7194908169548471, 1.2859924765981463))
add_parameter_sweep("sun", (0.9627690052874481, 1.6998823971640933))
add_parameter_sweep("m1p1", (1.2700977314900757, 2.2277696678863275))
add_parameter_sweep("m1p2", (1.5146762666232594, 2.6475155997806903))

def run_sweeps():
    for (sweep_dir, sweep_parameters) in sweeps.items():
        Parameter_Sweep(**sweep_parameters).generate_input_files(path("data", "parameter_sweeps", "synthesis", sweep_dir))
        run_multiplanet()

run_sweeps()