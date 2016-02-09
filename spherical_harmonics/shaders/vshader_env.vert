#version 400
uniform mat4 mvp_env;
in vec3 p_env;

out vec4 colorV;
void main () {
    gl_Position = mvp_env*vec4(p_env , 1.0);

    colorV = vec4(normalize(2+p_env),1.0);
}