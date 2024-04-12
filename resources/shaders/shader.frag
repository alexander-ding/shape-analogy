#version 330 core
out vec4 fragColor;

// Additional information for lighting
in vec4 normal_worldSpace;
in vec4 position_worldSpace;
in vec3 out_color;

uniform float red = 1.0;
uniform float green = 1.0;
uniform float blue = 1.0;
uniform float alpha = 1.0;

void main() {
    vec4 lightPos   = vec4(8.0, 8.0, -8.0 , 1.0);
    vec3 lightColor = vec3(1.0f, alpha, 0.0f);
    vec4 lightDir   = normalize(-lightPos + position_worldSpace);
    float c = 1; // clamp(dot(-normal_worldSpace, lightDir), 0, 1);

    fragColor = vec4(out_color[0] * c * lightColor[0], out_color[1] * c * lightColor[0], out_color[2] * c * lightColor[0], 1);
}
