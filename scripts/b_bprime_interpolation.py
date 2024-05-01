import numpy as np
import trimesh
import os
from scipy import interpolate
import numpy as np
from scipy.interpolate import interp1d

TIMESTEP = 40
MODEL_NAME = "" # your model name here (e.g. yoda)
B_PRIME = "" # path to b_prime
B_MESH = "" # path to b

# load in b and bprime
b_prime = trimesh.load(B_PRIME, maintain_order=True)
b = trimesh.load(B_MESH, maintain_order=True)
b_copy = b.copy()

# convert verts to numpy array 
b_prime_verts = np.array(b_prime.vertices)
b_verts = np.array(b.vertices)
curr_verts = b_verts

def delete_files_in_directory(directory_path):
   try:
     files = os.listdir(directory_path)
     for file in files:
       file_path = os.path.join(directory_path, file)
       if os.path.isfile(file_path):
         os.remove(file_path)
     print("All files deleted successfully.")
   except OSError:
     print("Error occurred while deleting files.")

folder = "" # folder to store each timestep's obj in 

# clear folder 
delete_files_in_directory(folder)

diff = (b_prime_verts - b_verts) / TIMESTEP

# interpolate
for i in range(0, TIMESTEP):
    b_copy.vertices += diff
    b_copy.export(f'{folder}/{MODEL_NAME}_{i}_.obj', include_texture = False)

