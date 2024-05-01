import trimesh
import os
from pathlib import Path
import numpy as np

TIMESTEP = 40
MODEL_NAME = ""  # your model name here (e.g. yoda)

input_dir = Path(__file__).parent.parent / 'interpolated'
output_dir = Path(__file__).parent.parent / 'output'
output_dir.mkdir(parents=True, exist_ok=True)

for f in output_dir.glob('*.obj'):
    os.remove(f)

steps = [
    trimesh.load(input_dir / 'sphere.obj'),
    trimesh.load(input_dir / 'd8.obj'),
    trimesh.load(input_dir / 'sphere.obj'),
    trimesh.load(input_dir / 'cylinder.obj'),
    trimesh.load(input_dir / 'sphere.obj'),
    trimesh.load(input_dir / 'tetrahedron.obj'),
    trimesh.load(input_dir / 'sphere.obj'),
    trimesh.load(input_dir / 'cylinder.obj'),
    trimesh.load(input_dir / 'sphere.obj'),
]

frame = 0
for i in range(1, len(steps)):
    after = steps[i]
    before = steps[i - 1]
    diff = (after.vertices - before.vertices) / TIMESTEP

    for i in range(TIMESTEP):
        before.export(output_dir / f'{frame}.obj')
        before.vertices += diff
        frame += 1

steps[-1].export(output_dir / f'{frame}.obj')
