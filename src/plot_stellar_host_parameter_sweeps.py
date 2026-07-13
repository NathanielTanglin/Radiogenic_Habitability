import nbformat
from paths import path

notebook_file = path("src", "plot_stellar_host_parameter_sweeps.ipynb")
nb = nbformat.read(notebook_file, as_version=nbformat.NO_CONVERT)

for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        exec(cell['source'])