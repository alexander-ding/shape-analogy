#version 330 core
out vec4 fragColor;

in vec4 normal_worldSpace;
in vec4 position_worldSpace;

uniform int   wire  = 0;
uniform float red   = 1.0;
uniform float green = 1.0;
uniform float blue  = 1.0;
uniform float alpha = 1.0;

void main() {
    // Do lighting in camera space
    vec4 lightPos   = vec4(0, 1000.0, 0 , 1.0);
    vec4 lightDir   = normalize(-lightPos + position_worldSpace);
    float c = clamp(dot(-normal_worldSpace, lightDir) + 0.2, 0, 1);

    fragColor = vec4(vec3(1, 1, 1) * c, 1);
}
