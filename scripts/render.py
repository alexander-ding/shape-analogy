import os
import numpy as np
import pyrender
import trimesh
import matplotlib.pyplot as plt
import shutil

model = "" # your model name here (e.g. "yoda")

def render_images(obj_filepath, ctr, path_to_pngs):
    model_trimesh = trimesh.load(obj_filepath)
    mesh = pyrender.Mesh.from_trimesh(model_trimesh)

    # compose scene
    scene = pyrender.Scene(ambient_light=[.1, .1, .3], bg_color=[0, 0, 0])
    camera = pyrender.PerspectiveCamera( yfov=np.pi / 3.0)
    light = pyrender.DirectionalLight(color=[1,1,1], intensity=2e3)

    scene.add(mesh, pose=  np.eye(4))
    scene.add(light, pose=  np.eye(4))

    # depending on the size of the mesh you load in, you'll need to change this
    # most likely just the translation column
    # +x points to the right, +y points up, +z points toward you (the screen)
    scene.add(camera, pose=[[ 1,  0,  0,  0],
                            [ 0,  1, 0, 0.5],
                            [ 0,  0, 1,  4],
                            [ 0,  0,  0,  1]])

    # render scene
    r = pyrender.OffscreenRenderer(512, 512)
    color, _ = r.render(scene)

    fig = plt.figure(figsize=(8,8), frameon= False)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    plt.imshow(color)
    filename = f'{path_to_pngs}/{ctr}.png'
    plt.savefig(filename)
    plt.close()
    # delete renderer 
    r.delete()


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
    

if __name__ == "__main__":
    path_to_objs = "" # your path to each timestep's obj files 
    path_to_pngs = "" # your path to each timestep's render 
    files = sorted(os.listdir(path_to_objs), key=lambda x: int(x.split('_')[1]))
    file_paths = files
    # clear existing images 
    delete_files_in_directory(path_to_pngs)
    
    ctr = 0
    for obj_file in file_paths: 
        try:
            render_images(obj_file, ctr, path_to_pngs)
            ctr += 1
        except Exception as e: # pyrender is kind of jank and throws a module not found error sometimes so just keep trying until it goes through lol
            render_images(obj_file, ctr, path_to_pngs)
            ctr += 1

    # this gives you renders going from sphere -> deformed 1
    # if you want to string together interpolations so it goes sphere -> deformed 1 -> sphere -> deformed 2 -> etc, 
    # first render images from sphere -> deformed 1, then copy those images in reverse order so you get sphere -> deformed 1 -> sphere images
    # that's what this code snippet does:

    # files = sorted(os.listdir(path_to_pngs), key=lambda x: int(x.split('.')[0]), reverse=True)
    # file_paths = files

    # for file in file_paths:
    #    source_file =  open(f'{path_to_pngs}/{file}', 'rb')
    #    destination_file = open(f'{path_to_pngs}/{ctr}.png', 'wb')
    #    shutil.copyfileobj(source_file, destination_file)
    #    ctr += 1

    # do this for all the deformed shapes you want 
    
    # have all the pics in one folder in the correct order, then string those together using ffmpeg:
    # ffmpeg -stream_loop 1 -framerate 48 -i "%d.png" -pix_fmt yuv420p output.mp4
    # changing # of loops and fps as needed 


