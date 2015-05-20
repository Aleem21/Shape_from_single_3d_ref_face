function [ depth ] = estimate_depth( N_ref_in, alb_ref, im, z_ref, sh_coeff,lambda )
%ESTIMATE_DEPTH Summary of this function goes here
%   Detailed explanation goes here
l0 = sh_coeff(1); l1 = sh_coeff(2); l2 = sh_coeff(3); l3 = sh_coeff(4);
[r_face,c_face] = find(~isnan(N_ref_in));

boundary = zeros(size(N_ref_in));
b = bwboundaries(~isnan(N_ref_in));b = b{1};
boundary(sub2ind(size(boundary),b(:,1),b(:,2))) = 1;
face_inds = sub2ind(size(im),r_face,c_face);
I = im(face_inds);
rho_ref = alb_ref(face_inds);
N_ref = N_ref_in(face_inds);
t1 = r_face+1;
t3 = c_face+1;

n_zs = numel(I);
A = zeros(n_zs, n_zs+1);
t1inds = sub2ind_my([n_zs,n_zs+1], 1:n_zs, t1);
t3inds = sub2ind_my([n_zs,n_zs+1], 1:n_zs, t3);

A(t1inds) = l1;
A(t3inds) = l2;
try
for i=1:n_zs
    
    A(i,i) = -l1 -l2;
    %find x+1,y
    if r_face(i)+1>size(im,1)
        elem_num = find_elements([r_face c_face],r_face(i),c_face(i));
        A(i, elem_num ) = l1;
    elseif ~(boundary(r_face(i)+1,c_face(i)))
        elem_num = find_elements([r_face c_face],r_face(i)+1,c_face(i));
        A(i, elem_num ) = l1;
    else
        A(i,end) = A(i,end) + l1*z_ref(r_face(i)+1,c_face(i));
    end
    if c_face(i)+1>size(im,2)
        elem_num = find_elements([r_face c_face],r_face(i),c_face(i));
        A(i, elem_num) = l2;
    elseif ~(boundary(r_face(i),c_face(i)+1))
        elem_num = find_elements([r_face c_face],r_face(i),c_face(i)+1);
        A(i, elem_num) = l2;
    else
        A(i,end) = A(i,end) + l2*z_ref(r_face(i),c_face(i)+1);
    end
    
end
catch err
    disp(err.message)
end
rhs = I - rho_ref*l0 + rho_ref./N_ref*l3;
depth = A\rhs;
