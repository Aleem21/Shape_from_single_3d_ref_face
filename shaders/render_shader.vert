#version 400
in vec3 p;
in vec3 rgb;
out vec4 rgb_f;
void main () {
    gl_Position = vec4(p, 1.0);
    rgb_f = vec4(rgb,1);
}