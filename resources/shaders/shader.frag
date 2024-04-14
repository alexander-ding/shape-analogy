#version 330 core
out vec4 fragColor;

// Additional information for lighting
in vec4 normal_worldSpace;
in vec4 position_worldSpace;
in vec3 out_color;

void main() {
    vec4 lightPos   = vec4(8.0, 8.0, -8.0 , 1.0);
    vec4 lightDir   = normalize(-lightPos + position_worldSpace);
    float c = 1; // clamp(dot(-normal_worldSpace, lightDir), 0, 1);

    fragColor = vec4(out_color * c, 1);
}
