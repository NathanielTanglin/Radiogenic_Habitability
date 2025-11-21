import nbformat

# Note: add modified multiplanet binary into project directory and modify the PATH variable.

# Runs the Jupyter notebook with all the parameter sweep scripts.
nb = nbformat.read('run_parameter_sweeps.ipynb', as_version=nbformat.NO_CONVERT)

for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        exec(cell['source'])