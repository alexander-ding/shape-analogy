import os
import numpy as np
import pyrender
import trimesh
import matplotlib.pyplot as plt
from pathlib import Path


def look_at(position, target, up_vector=np.array([0, 1, 0])):
    # Calculate forward, right, and up vectors
    direction = position - target
    direction = direction / np.linalg.norm(direction)

    right_vector = np.cross(up_vector, direction)
    right_vector = right_vector / np.linalg.norm(right_vector)

    new_up_vector = np.cross(direction, right_vector)
    new_up_vector = new_up_vector / np.linalg.norm(new_up_vector)

    # if np.dot(new_up_vector, up_vector) < 0:
    #     new_up_vector = -new_up_vector

    # Construct the rotation matrix
    rotation_matrix = np.column_stack((
        right_vector,
        new_up_vector,
        direction
    ))

    # Construct the full transformation matrix
    pose = np.eye(4)
    pose[:3, :3] = rotation_matrix
    pose[:3, 3] = position

    return pose


def render_image(obj_path, output_path):
    obj = trimesh.load(obj_path)
    obj.apply_transform(trimesh.transformations.rotation_matrix(
        np.pi, [0, 1, 0], [0, 0, 0]))
    obj.vertices -= obj.vertices.mean(axis=0)
    mesh = pyrender.Mesh.from_trimesh(obj)

    # compose scene
    scene = pyrender.Scene(ambient_light=[.1, .1, .3], bg_color=[0, 0, 0])
    # camera = pyrender.OrthographicCamera()
    camera = pyrender.PerspectiveCamera(yfov=np.pi / 3.0)
    light = pyrender.DirectionalLight(color=[1, 1, 1], intensity=2e3)

    scene.add(mesh, pose=np.eye(4))
    scene.add(light, pose=np.eye(4))

    pos = np.array([1.5, 1, 3])
    look = np.array([0, 0, 0])
    pose = look_at(pos, look)
    # Set up the scene
    scene.add(camera, pose=pose)

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
