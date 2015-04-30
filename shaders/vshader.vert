#version 400
in vec3 p;
out vec4 colorV;
void main () {
    gl_Position = vec4(p, 1.0);
    colorV = vec4(1,1,p.z,1);
}