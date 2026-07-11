import nbformat
import os
from paths import path

notebook_file = path("src", "run_GKM_host_parameter_sweeps.ipynb")

# Runs the Jupyter notebook with the parameter sweep code for the GV-, KV-, and MV-systems.
nb = nbformat.read(notebook_file, as_version=nbformat.NO_CONVERT)

for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        exec(cell['source'])