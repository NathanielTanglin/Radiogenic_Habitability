## Description

This repository stores the code used in Tanglin & Becker (2026).

## Instructions

WINDOWS NOT SUPPORTED. I recommend using Windows Subsystem for Linux (WSL) if you have a Windows OS.

**To install the code, do the following:**

1. Clone this repository (`git clone https://github.com/NathanielTanglin/Radiogenic_Habitability.git`). The project directory must be named Radiogenic_Habitability. If you change the name, you must modify `paths.py`.
2. In a terminal, run `python setup.py develop` in the `modded_VPLanet` directory. This will install VPLanet onto your machine with our modifications to `atmesc.c` included. (See the **Tips** section below to hopefully make this easier on yourself.)
3. You will need the `nbformat` python package. Install this using whatever package manager you use (e.g. `pip install nbformat` or `conda install nbformat`). This module is necessary for running `run_GKM_host_parameter_sweeps.py`, `run_all_parameter_sweeps.py`, and `plot_stellar_host_parameter_sweeps.py`.

**To run the simulations (i.e. populate the data directory), you may do any of the following:**

| Action | Description |
| ----------- | ----------- |
| `python run_GKM_host_parameter_sweeps.py` | Runs the simulations corresponding to Figure 4 of the paper. These simulations take the longest, and WILL TAKE UP TO 90 GB OF SPACE. |
| `python run_stellar_host_mass_parameter_sweeps.py` | Runs the simulations corresponding to Figures 5 and 6. These simulations will finish in less than 20 minutes given 8 cores on a reasonable computer, and will take about 4 GB of space. |
| `python run_all_parameter_sweeps.py` | Runs _fingers crossed_ all the above simulations. |

**To visualize data (i.e. populate the plots directory), you may do any of the following:**

| Action | Description |
| ----------- | ----------- |
| `python plot_planet_B_F_AP.py` | Plots Figure 1 |
| `python plot_thermal_evolution.py` | Plots Figure 2 |
| `python plot_abundances.py` | Plots Figure 3 |
| `python plot_GKM_host_parameter_sweeps.py` | Plots Figures 4a-f |
| `python plot_stellar_host_parameter_sweeps.py` | Plots Figures 5a-f and Figure 6 |

## Tips

**I would recommend using** [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) to create a virtual environment with all the VPLanet dependencies installed via `modded_VPLanet/environment.yml`:

```
cd modded_VPLanet
conda create -n vplanet_env --environment-specifier environment.yml
conda activate vplanet_env
python setup.py develop
conda install nbformat
```

**If you want to run this code on a computing cluster** in which you submit jobs, you will need to edit `run_GKM_host_parameter_sweeps.ipynb` by setting `num_cores="[your number of CPUs]"` (make sure the number is in quotes!) to match whatever you have on your job script. You will need to do the same for `run_stellar_host_masss_parameter_sweeps.py`.
