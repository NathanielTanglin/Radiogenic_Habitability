import os

project_dir = 'Radiogenic_Habitability_Replicable'
root = os.getcwd().split(project_dir)[0] + project_dir

def path(*args):
    return os.path.join(root, os.sep.join(args))