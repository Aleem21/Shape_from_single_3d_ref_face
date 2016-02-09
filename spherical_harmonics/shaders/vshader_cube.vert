#version 400
uniform mat4 mvp_env;
uniform mat3 rot_env;
in vec3 p_env;

out vec3 TexCoords;
void main() {
    
    gl_Position = mvp_env*vec4(p_env, 1.0);
    TexCoords =  rot_env*p_env;
}