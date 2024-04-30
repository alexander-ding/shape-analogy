from argparse import ArgumentParser
from pathlib import Path
import trimesh
import math
import numpy as np
import os

root_dir = Path(__file__).parent.parent
output_dir = root_dir / 'output'
for o in output_dir.glob('*.obj'):
    os.remove(o)
parser = ArgumentParser("Interpolate between two meshes")
parser.add_argument("mesh", type=Path)
parser.add_argument("-n", type=int, default=60)
parser.add_argument("-o", type=Path, default=output_dir)
args = parser.parse_args()
mesh = trimesh.load(args.mesh)

bounds = mesh.bounds
desired_scale = 2 * math.sqrt(3)
intersection_scale = 2 / 4

# fit mesh within [-0.5, 0.5]
scale_factor = intersection_scale / max(bounds[1] - bounds[0])
mesh.apply_scale(scale_factor)
translation = -mesh.bounds.mean(axis=0)
mesh.apply_translation(translation)

output_dir = args.o
output_dir.mkdir(exist_ok=True, parents=True)
sphere_mesh = trimesh.creation.icosphere(subdivisions=6)

vertices = sphere_mesh.vertices
dirs = -trimesh.util.unitize(vertices)
locations, index_ray, index_tri = mesh.ray.intersects_location(
    ray_origins=vertices,
    ray_directions=dirs
)

# Calculate distances and find the closest intersection for each ray
target_vertices = np.full_like(vertices, np.nan)  # Initialize with NaNs

for i in range(vertices.shape[0]):
    mask = index_ray == i
    if np.any(mask):
        hits = locations[mask]
        distances = np.linalg.norm(hits - vertices[i], axis=1)
        closest_hit_index = np.argmin(distances)
        closest_point = hits[closest_hit_index]
        closest_tri = mesh.faces[index_tri[closest_hit_index]]

        # Store the point
        target_vertices[i] = closest_point
target_vertices = target_vertices / intersection_scale
original_vertices = sphere_mesh.vertices
delta = (target_vertices - vertices) / args.n
for i in range(args.n + 1):
    sphere_mesh.vertices = original_vertices + delta * i
    sphere_mesh.export(output_dir / f"{i}.obj")
