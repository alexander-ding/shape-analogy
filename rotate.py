import math

def rotate_vertex(x, y, z, angle):
    """Rotate the vertex around the Y-axis by the given angle in degrees."""
    rad = math.radians(angle)
    # Rotate around Y-axis
    x_new = z * math.sin(rad) + x * math.cos(rad)
    z_new = z * math.cos(rad) - x * math.sin(rad)
    return x_new, y, z_new

def rotate_obj_file(input_path, output_path, angle):
    with open(input_path, 'r') as file:
        lines = file.readlines()

    with open(output_path, 'w') as file:
        for line in lines:
            if line.startswith('v '):
                parts = line.split()
                x, y, z = map(float, parts[1:4])
                x, y, z = rotate_vertex(x, y, z, angle)
                file.write(f'v {x} {y} {z}\n')
            else:
                file.write(line)

# Example usage
input_obj_path = 'meshes/moai.obj'  # Path to your original OBJ file
output_obj_path = 'meshes/moai-rotated.obj'  # Path to save the rotated OBJ file
rotate_angle = 180  # Angle in degrees

rotate_obj_file(input_obj_path, output_obj_path, rotate_angle)
