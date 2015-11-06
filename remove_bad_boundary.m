function [ face ] = remove_bad_boundary( face )
%REMOVE_BAD_BOUNDARY Remove such points from face on thich boundary
%conditions would not be possible to compute eventually, i.e. no face
%pixels on (top AND bottom) OR (left AND right)
face = double(face);
hor = [1 0 1];
vert = hor';

f_vert = (conv2(face,vert,'same').*face)>0;
f_hor = (conv2(face,hor,'same').*face)>0;
if sum(sum(f_vert))~=sum(face(:)) || sum(sum(f_hor))~=sum(face(:))
    face = f_vert.*f_hor;
    face = remove_bad_boundary(face);
else
    face = f_vert.*f_hor;
end
face = face>0;
end

