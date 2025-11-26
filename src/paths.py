import os

# Gets the current working directory, splits it by the file separator,
# and takes the second to last entry (the project directory).
project_dir = os.getcwd().split(os.sep)[-2]
root = os.getcwd().split(project_dir)[0] + project_dir

def path(*args):
    return os.path.join(root, os.sep.join(args))