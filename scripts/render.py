import os
import numpy as np
import pyrender
import trimesh
import matplotlib.pyplot as plt
from pathlib import Path


def render_image(obj_path, output_path):
    obj = trimesh.load(obj_path)
    obj.apply_transform(trimesh.transformations.rotation_matrix(
        np.pi, [0, 1, 0], [0, 0, 0]))
    obj.vertices -= obj.vertices.mean(axis=0)
    mesh = pyrender.Mesh.from_trimesh(obj)

    # compose scene
    scene = pyrender.Scene(ambient_light=[.1, .1, .3], bg_color=[0, 0, 0])
    camera = pyrender.PerspectiveCamera(yfov=np.pi / 3.0)
    light = pyrender.DirectionalLight(color=[1, 1, 1], intensity=2e3)

    scene.add(mesh, pose=np.eye(4))
    scene.add(light, pose=np.eye(4))

    # depending on the size of the mesh you load in, you'll need to change this
    # most likely just the translation column
    # +x points to the right, +y points up, +z points toward you (the screen)
    scene.add(camera, pose=[[1,  0,  0,  0],
                            [0,  1, 0, 0.5],
                            [0,  0, 1,  4],
                            [0,  0,  0,  1]])

    # render scene
    r = pyrender.OffscreenRenderer(512, 512)
    color, _ = r.render(scene)

    fig = plt.figure(figsize=(8, 8), frameon=False)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    plt.imshow(color)
    plt.savefig(output_path)
    plt.close()
    # delete renderer
    r.delete()


if __name__ == "__main__":
    output_dir = Path(__file__).parent.parent / 'render'
    input_dir = Path(__file__).parent.parent / 'output'
    all_files = sorted(input_dir.glob("*.obj"),
                       key=lambda x: int(x.name.split('.')[0]))

    for obj_file in all_files:
        frame = obj_file.name.split('.')[0]
        print(obj_file)
        while True:
            try:
                render_image(obj_file, output_dir /
                             f"{frame}.png")
                break
            except Exception as e:
                print(e)
                continue

    # have all the pics in one folder in the correct order, then string those together using ffmpeg:
    # ffmpeg -stream_loop 1 -framerate 48 -i "%d.png" -pix_fmt yuv420p output.mp4
    # changing # of loops and fps as needed
