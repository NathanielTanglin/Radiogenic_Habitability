import os

# Finds the absolute path to the root project directory.
root = os.path.join(os.getcwd().split("Radiogenic_Habitability")[0], "Radiogenic_Habitability") #os.getcwd().split("Radiogenic_Habitability")[0]

def path(*args):
    return os.path.join(root, os.sep.join(args))