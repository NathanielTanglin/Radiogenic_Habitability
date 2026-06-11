**Instructions:**
1. Go to the modded_VPLanet directory and run "python setup.py develop" in a Linux terminal (you can use Windows Subsystem for Linux on a Windows computer). This will install VPLanet onto your machine.
2. Then run "pip install nbformat", or whatever package manager you use instead of pip. The nbformat module is necessary for calling the run_all_parameter_sweeps.py script.
3. Finally, call "python run_all_parameter_sweeps.py" in terminal. This will *fingers crossed* run all the relevant simulations and populate the data directory with results.
4. To visualize the data, call any of plot_thermal_evolution.py, plot_synthesis.py, plot_contours.py, or plot_abundances.py. This should produce all the relevant figures.
