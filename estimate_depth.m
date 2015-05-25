function [ depth ] = estimate_depth( N_ref_in, alb_ref, im, z_ref, sh_coeff,lambda1 )
%ESTIMATE_DEPTH Summary of this function goes here
%   Detailed explanation goes here
l0 = sh_coeff(1); l1 = sh_coeff(2); l2 = sh_coeff(3); l3 = sh_coeff(4);
n_rows = size(im,1);
n_cols = size(im,2);
[b_out,b_in] = find_boundary(~isnan(N_ref_in),false);
[b_out_full,b_in_full] = find_boundary(~isnan(N_ref_in),true);
face = ~isnan(N_ref_in);
in_face = ~isnan(N_ref_in)-b_out_full;
[r_face,c_face] = find(face);

%contour normals
[ncx,ncy] = find_countour_normal(b_out);
% [b_out,b_in] = find_boundary(~isnan(N_ref_in),0);
ncx(~b_out) = NaN;
ncy(~b_out) = NaN;
face_inds = sub2ind(size(im),r_face,c_face);
ncx = ncx(face_inds);
ncy = ncy(face_inds);
I = im(face_inds);
rho_ref = alb_ref(face_inds);
N_ref = N_ref_in(face_inds);
t1 = r_face+1;
t3 = c_face+1;

n_zs = numel(I);
% A = zeros(n_zs, n_zs+1);
tic
n_boundary = sum(b_out(:));
A = sparse(n_zs+n_boundary, n_zs);

% t1inds = sub2ind_my([n_zs,n_zs+1], 1:n_zs, t1);
% t3inds = sub2ind_my([n_zs,n_zs+1], 1:n_zs, t3);
%s
% A(t1inds) = l1;
% A(t3inds) = l2;
counter = n_zs+1;
f = figure;
bad_boundary = [];
[fr,fc] = meshgrid(1:size(face,1),1:size(face,2));
inds = find_elements([r_face c_face],fr',fc');
try
for i=1:n_zs
    % if its boundary use boundary conditions
    if b_out(r_face(i),c_face(i))
        % mark it to remove later
        A(i,1) = NaN;
        %-- for x (rows)
        %if location below is out of bound or not inside face,
        %subtract last entry from current
        if (r_face(i)+1)>n_rows || ~face(r_face(i)+1,c_face(i))
%             neg = find_elements([r_face c_face],r_face(i)-1,c_face(i));
            neg = inds(r_face(i)-1,c_face(i));
            pos = i;
            % else subtract current entry from next one
        else
            neg = i;
%             pos = find_elements([r_face c_face],r_face(i)+1,c_face(i));
            pos = inds(r_face(i)+1,c_face(i));
        end
        % remove boundary points on which boundary conditions can not be
        % computed
        if sum([neg pos]==0)
            bad_boundary(end+1) = i;
        else
            A(counter,neg) = A(counter,neg)-ncx(i);
            A(counter,pos) = A(counter,pos)+ncx(i);
        end
        % -- for y (cols)
        % if location on right is out of bound or not inside face,
        % subtract last entry from current
        if (c_face(i)+1)>n_cols || ~face(r_face(i),c_face(i)+1)
%             neg = find_elements([r_face c_face],r_face(i),c_face(i)-1);
            neg = inds(r_face(i),c_face(i)-1);
            pos = i;
            % else subtract current entry from next one
        else
            neg = i;
%             pos = find_elements([r_face c_face],r_face(i),c_face(i)+1);
            pos = inds(r_face(i),c_face(i)+1);
        end
        if sum([neg pos]==0)
            bad_boundary(end+1) = i;
        else
            A(counter,neg) = A(counter,neg)-ncy(i);
            A(counter,pos) = A(counter,pos)+ncy(i);
        end
        counter = counter + 1;
    else
        %if b_in(r_face(i),c_face(i))
        %	neg = find_elements([r_face c_face],r_face(i)-1,c_face(i));
        %	pos = find_elements([r_face c_face],r_face(i)+1,c_face(i));
        %	A(counter,neg) = A(counter,neg)-ncx(i);
        %	A(counter,pos) = A(counter,pos)+ncx(i);
        %	neg = find_elements([r_face c_face],r_face(i)-1,c_face(i));
        %	pos = find_elements([r_face c_face],r_face(i)+1,c_face(i));
        %	A(counter,neg) = A(counter,neg)-ncy(i);
        %	A(counter,pos) = A(counter,pos)+ncy(i);
        %	counter = counter + 1;
        %	% if its an inside the face pixel, use data term constraint
        %end
        A(i,i) = -l1 -l2;
        
%         elems = find_elements([r_face c_face],[1 0]+r_face(i),[0 1]+c_face(i));
        elems = inds(sub2ind(size(face),[1 0]+r_face(i),[0 1]+c_face(i)));
        A(i, elems ) = [l1 l2];
    end
end
catch err
    disp 1
end

% regularization term
reg = 1;
if (reg)
[r_inface,c_inface] = find(in_face);
inface_inds = sub2ind(size(im),r_inface,c_inface);
n_in_zs = numel(inface_inds);
sz = 3; dev = 2;
gauss = fspecial('gaussian',sz,dev);
% gauss = [-1 -1  -1; -1 8 -1; -1 -1 -1];
rhs_reg = lambda1*(z_ref - conv2(z_ref,gauss,'same'));
% rhs_reg = lambda1*(conv2(z_ref,gauss,'same'));
rhs_reg = rhs_reg(inface_inds);
[boxc, boxr] = meshgrid(-1:1,-1:1);

gauss = -gauss;
gauss(2,2) = 1+gauss(2,2);
gaussVec = lambda1*gauss(:);
% gaussVec = -lambda1*gauss(:);
try
for i=1:n_in_zs
     elems3x3 = inds(sub2ind(size(face),boxr(:)+r_inface(i),boxc(:)+c_inface(i)));
%     elems3x3 = find_elements([r_face c_face],boxr(:)+r_inface(i),boxc+c_inface(i));
    A(end+1,elems3x3) = gaussVec;
%     mid_i = elems3x3(5);
%     A(end,mid_i) = A(end,mid_i) + 1;
%     A(end,:) = A(end,:)*lambda1;
end
catch err
    disp 1
end
end
% for i=1:n_zs
%     % find elements in order NW, N, NE, W, MID, E, SW, S, SE
%     inds = find_elements([r_face c_face],[-1 -1 -1 0 0 0 1 1 1]+r_inface(i),[-1 0 1 -1 0 1 -1 0 1]+c_inface(i));
% end
toc
rhs = (I./rho_ref - l0).*N_ref + l3;
rhs(end+1:end+n_boundary) = 0;
rhs = [rhs; rhs_reg];
% remove rows in A where data term was not computed. Also remove
% corresponsing image observations (i.e. LHS of data term)
bad_rows = isnan(A(:,1));
A(bad_rows,:) = [];
rhs(bad_rows) = [];

% remove columns of A for which depth can not be determined because that
% depth has not been constrained by any boudnary condition as the gradient
% on that pixel could not be computed
A(:,bad_boundary) = [];
face_inds(bad_boundary) = [];


z = A\rhs;
depth = zeros(size(N_ref_in));
depth(face_inds) = z;
offset = depth(inface_inds(1)) - z_ref(inface_inds(1));
depth = depth-offset;
depth(~face) = NaN;
figure; surf(depth,'edgealpha',0);
z2 = z_ref(face_inds);
err2 = abs(A*z2-rhs);
err = abs(A*z-rhs);
depth2 = depth;
depth2(face_inds) = err(1:numel(face_inds));
figure;imagesc(depth2)
