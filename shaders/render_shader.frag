#version 400
out vec4 frag_colour;
in vec4 rgb_f;
void main () {
   frag_colour  = rgb_f;
}