function [ depth ] = estimate_depth( N_ref_in, alb_ref, im, z_ref, sh_coeff,lambda1,reg_type,z_gnd,talk)
%ESTIMATE_DEPTH Summary of this function goes here
%   Detailed explanation goes here
if nargin < 8
    talk = 1;
end
a0 = pi;
a1 = 2*pi/sqrt(3);
a2= 2*pi/sqrt(8);

l0 = sh_coeff(1); l1 = sh_coeff(2); l2 = sh_coeff(3); l3 = sh_coeff(4);
n_rows = size(im,1);
n_cols = size(im,2);

face = ~isnan(N_ref_in);
face = remove_bad_boundary(face);
[b_out,b_in] = find_boundary(face,false);
[b_out_full,b_in_full] = find_boundary(face,true);
in_face = face-b_out_full;
[r_face,c_face] = find(face);

%contour normals
[ncx,ncy] = find_countour_normal(b_out);
% [b_out,b_in] = find_boundary(~isnan(N_ref_in),0);
ncx(~b_out) = NaN;
ncy(~b_out) = NaN;
face_inds = sub2ind(size(im),r_face,c_face);

ncx = -ncx(face_inds);
ncy = -ncy(face_inds);
I = im(face_inds);
rho_ref = alb_ref(face_inds);
N_ref = N_ref_in(face_inds);
% t1 = r_face+1;
% t3 = c_face+1;

n_zs = numel(I);
% A = zeros(n_zs, n_zs+1);
if talk
    tic
end
n_boundary = sum(b_out(:));
% A = sparse(n_zs+n_boundary, n_zs);
A = sparse([],[],[],n_zs+n_boundary, n_zs,sum(in_face(:))*3 + sum(b_out(:))*4 + sum(in_face(:))*9);

% t1inds = sub2ind_my([n_zs,n_zs+1], 1:n_zs, t1);
% t3inds = sub2ind_my([n_zs,n_zs+1], 1:n_zs, t3);
%s
% A(t1inds) = l1;
% A(t3inds) = l2;
counter = n_zs+1;
% f = figure;
bad_boundary = [];
[fr,fc] = meshgrid(1:size(face,1),1:size(face,2));
inds = find_elements([r_face c_face],fr',fc');
bad_bnd_cur = 0;
dark_factor = 1;
try
for i=1:n_zs
    % if its boundary use boundary conditions
    
    if b_out(r_face(i),c_face(i))
        % mark it to remove later
        A(i,1) = NaN;
        %-- for x (rows)
        %if location below is out of bound or not inside face,
        %subtract last entry from current
        if (r_face(i)-1)<1 || ~face(r_face(i)-1,c_face(i))
%             neg = find_elements([r_face c_face],r_face(i)-1,c_face(i));
            neg = inds(r_face(i)+1,c_face(i));
            pos = i;
            % else subtract current entry from next one
        else
            neg = i;
%             pos = find_elements([r_face c_face],r_face(i)+1,c_face(i));
            pos = inds(r_face(i)-1,c_face(i));
        end
        % remove boundary points on which boundary conditions can not be
        % computed
        if sum([neg pos]==0)
            bad_boundary(end+1) = i;
            bad_bnd_cur = 1;
        else
            A(counter,neg) = A(counter,neg)-ncx(i);
            A(counter,pos) = A(counter,pos)+ncx(i);
        end
        % -- for y (cols)
        % if location on right is out of bound or not inside face,
        % subtract last entry from current
        if (c_face(i)-1)<1 || ~face(r_face(i),c_face(i)-1)
%             neg = find_elements([r_face c_face],r_face(i),c_face(i)-1);
            neg = inds(r_face(i),c_face(i)+1);
            pos = i;
            % else subtract current entry from next one
        else
            neg = i;
%             pos = find_elements([r_face c_face],r_face(i),c_face(i)+1);
            pos = inds(r_face(i),c_face(i)-1);
        end
        if sum([neg pos]==0)
            bad_boundary(end+1) = i;
            bad_bnd_cur = 1;
        else
            A(counter,neg) = A(counter,neg)-ncy(i);
            A(counter,pos) = A(counter,pos)+ncy(i);
        end
        if bad_bnd_cur
            A(counter,:) = 0;
            bad_bnd_cur = 0;
        else
            counter = counter + 1;
        end
    else
        factor = alb_ref(r_face(i),c_face(i))/N_ref_in(r_face(i),c_face(i));
        if I(i)<0.01
            factor = factor/dark_factor;
        end
        %         const = alb_ref(r_face(i),c_face(i))*l0 - factor*l3;
%         elems = inds(sub2ind(size(face),[0 +1]+r_face(i),[+1 0]+c_face(i)));
        A(i,i) = (-l1 -l2)*factor;
%         A(i,elems) = [-l1 -l2]*factor;
%                 elems = inds(sub2ind(size(face),[-1 0]+r_face(i),[0 -1]+c_face(i)));
        %         A(i, elems ) = -factor*[l1 l2];
        %         elems = find_elements([r_face c_face],[1 0]+r_face(i),[0 1]+c_face(i));
        elems = inds(sub2ind(size(face),[0 -1]+r_face(i),[-1 0]+c_face(i)));
        A(i, elems ) = factor*[l1 l2];
        %         A(i,end) = const;  

    end
end
catch err
    error('some error')
end
A(end-numel(bad_boundary)+1:end,:)=[];
% regularization term
% in_face = in_face-b_in;
[r_inface,c_inface] = find(in_face);

inface_inds = sub2ind(size(im),r_inface,c_inface);


reg = 1;
if (reg)
n_in_zs = numel(inface_inds);
sz =3; dev = 2;
gauss = fspecial('gaussian',sz,dev);
% gauss = [-1 -1  -1; -1 8 -1; -1 -1 -1];
% gauss = [0 -1  0; -1 4 0; 0 -1 0];
rhs_reg_mat = lambda1*(z_ref - conv2(z_ref,gauss,'same'));

% del me dz only
if strcmp(reg_type,'dz')
    rhs_reg_mat = lambda1*(z_ref);
end
% end del me dz only


% rhs_reg_mat = lambda1*(conv2(z_ref,gauss,'same'));
rhs_reg = rhs_reg_mat(inface_inds);
f_w = floor(sz/2);
[boxc, boxr] = meshgrid(-f_w:f_w,-f_w:f_w);

gauss = -gauss;
gauss(2,2) = 1+gauss(2,2);

%del me dz only
if strcmp(reg_type,'dz')
    gauss = gauss*0;
    gauss(2,2) = 1;
end
%end del me dz only

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
    error( 'error!!!')
end
end
% for i=1:n_zs
%     % find elements in order NW, N, NE, W, MID, E, SW, S, SE
%     inds = find_elements([r_face c_face],[-1 -1 -1 0 0 0 1 1 1]+r_inface(i),[-1 0 1 -1 0 1 -1 0 1]+c_inface(i));
% end
if talk
    toc
end
% rhs = (I./rho_ref - l0).*N_ref + l3;
rhs = (I - rho_ref*l0 + rho_ref./N_ref*l3);
rhs(I<0.01) = rhs(I<0.01)/dark_factor;
% rhs(rhs_is_boundary>0) = 0;
rhs(end+1:end+n_boundary-numel(bad_boundary)) = 0;
if reg
    rhs = [rhs; rhs_reg];
end
% remove rows in A where data term was not computed. Also remove
% corresponsing image observations (i.e. LHS of data term)

bad_rows = isnan(A(:,1));
A(bad_rows,:) = [];
rhs(bad_rows) = [];

constraint_inds = 1:n_zs;
constraint_inds(bad_rows) = [];
constraint_inds = sub2ind(size(face),r_face(constraint_inds),c_face(constraint_inds));

% remove columns of A for which depth can not be determined because that
% depth has not been constrained by any boudnary condition as the gradient
% on that pixel could not be computed
A(:,bad_boundary) = [];
face_inds(bad_boundary) = [];


z = A\rhs;
% zend = z(end);
% z = z(1:end-1)/z(end);
depth = NaN(size(N_ref_in));
depth(face_inds) = z;
offset = mean(depth(inface_inds)) - mean(z_ref(inface_inds));
depth = depth-offset;
depth(~face) = NaN;
% % figure; surf(depth,'edgealpha',0);

if nargin>7
    z2 = z_gnd(face_inds);
else
    z2 = z_ref(face_inds);
end
err_ref = (A*z2-rhs).^2;
err_est = (A*z-rhs).^2;
% depth2 = depth;
% depth2(face_inds) = err_est(1:numel(face_inds));
% depth22 = depth;
% depth22(face_inds) = err_ref(1:numel(face_inds));
s = 1; e = n_zs-n_boundary;
bnd_s = e+1;
bnd_e = e+n_boundary-numel(bad_boundary);
reg_s = bnd_e+1;
reg_e = bnd_e+n_in_zs;
if talk
    figure;plot(s:e+1,err_ref(s:e+1),bnd_s:bnd_e+1,err_ref(bnd_s:bnd_e+1),reg_s:reg_e,err_ref(reg_s:reg_e));title('error in ref depth')
    fprintf('Ground truth error = %d\n',sum(err_ref.^2));
    figure;plot(s:e+1,err_est(s:e+1),bnd_s:bnd_e+1,err_est(bnd_s:bnd_e+1),reg_s:reg_e,err_est(reg_s:reg_e));title('error in computed depth')
    fprintf('Estimated error = %d\n',sum(err_est.^2));
    face_cost_ref = NaN(size(face));
    face_cost_ref(constraint_inds) = err_ref(s:e);
    figure;imagesc(face_cost_ref);title('constraint cost, ref face')
    face_cost_est = NaN(size(face));
    face_cost_est (constraint_inds) = err_est(s:e);
    figure;imagesc(face_cost_est);title('constraint cost, ref est')
end
% figure;imagesc(depth);
% figure;imagesc(depth2)
