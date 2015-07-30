function [ xp,xn,yp,yn,pref] = find_grad_inds_boundary( r,c,face,sub2ind_face,pref_in )
%FIND_GRAD_INDS_BOUNDARY Summary of this function goes here
%   Detailed explanation goes here
[szr,szc] = size(face);
if pref_in(2)<0
    if (r+1)>szr || ~face(r+1,c)
        neg = sub2ind_face(r,c);          %me
        pos = sub2ind_face(r-1,c);        %up
        pref(2) = 1;
    else
        %if up is not accesable
        neg = sub2ind_face(r+1,c);          %down
        pos = sub2ind_face(r,c);          %me
        % else subtract current entry from next one
        pref(2) = -1;
    end
else
    if (r-1)<1 || ~face(r-1,c)
        %if up is not accesable
        neg = sub2ind_face(r+1,c);          %down
        pos = sub2ind_face(r,c);          %me
        % else subtract current entry from next one
        pref(2) = -1;
    else
        neg = sub2ind_face(r,c);          %me
        pos = sub2ind_face(r-1,c);        %up
        pref(2) = 1;
    end
end
yp = pos;
yn = neg;

if pref_in(1)<0
    if (c+1)>szc || ~face(r,c+1)
        neg = sub2ind_face(r,c);          %left
        pos = sub2ind_face(r,c-1);        %up
        pref(1) = 1;
    else
        neg = sub2ind_face(r,c+1);          %right
        pos = sub2ind_face(r,c);          %me
        pref(1) = -1;
    end
else
    if (c-1)<1 || ~face(r,c-1)
        %if up is not accesable
        neg = sub2ind_face(r,c+1);          %right
        pos = sub2ind_face(r,c);          %me
        pref(1) = -1;
        % else subtract current entry from next one
    else
        neg = sub2ind_face(r,c);          %me
        pos = sub2ind_face(r,c-1);        %left
        pref(1) = 1;
    end
end
xp = pos;
xn = neg;
end

