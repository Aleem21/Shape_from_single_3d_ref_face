function [ depth ] = estimate_depth_nonlin( N_ref_in, alb_ref, im, z_ref, sh_coeff,lambda1,reg_type,z_gnd,talk)
%ESTIMATE_DEPTH Summary of this function goes here
%   Detailed explanation goes here
if nargin < 8
    talk = 1;
end

l0 = sh_coeff(1); l1 = sh_coeff(2); l2 = sh_coeff(3); l3 = sh_coeff(4);

face = ~isnan(N_ref_in);
face = remove_bad_boundary(face)>0;
[b_out,b_in] = find_boundary(face,false);
[b_out_full,b_in_full] = find_boundary(face,true);
in_face = (face-b_out_full)>0;
[r_inface,c_inface] = find(in_face);
inface_inds = sub2ind(size(im),r_inface,c_inface);
[r_face,c_face] = find(face);

%contour normals
[ncy,ncx] = find_countour_normal(b_out_full);
% [b_out,b_in] = find_boundary(~isnan(N_ref_in),0);
ncx(~b_out_full) = NaN;
ncy(~b_out_full) = NaN;
face_inds = sub2ind(size(im),r_face,c_face);
face_inds2d = nan(size(face));
face_inds2d(face) = face_inds;
ncx = -0.0001*ncx(b_out_full);
ncy = -0.0001*ncy(b_out_full);
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

xp = inds(sub2ind(size(im),r_inface,c_inface-1));
xn = inds(sub2ind(size(im),r_inface,c_inface));

yp = inds(sub2ind(size(im),r_inface-1,c_inface));
yn = inds(sub2ind(size(im),r_inface,c_inface));
[r_bound,c_bound] = find(b_out_full);
yp_bound = [];
yn_bound = [];

xp_bound = [];
xn_bound = [];
for i=1:sum(b_out_full(:))
    
    if (r_bound(i)-1)<1 || ~face(r_bound(i)-1,c_bound(i))
        %if up is not accesable
        neg = inds(r_bound(i)+1,c_bound(i));          %down
        pos = inds(r_bound(i),c_bound(i));          %me
        % else subtract current entry from next one
    else
        neg = inds(r_bound(i),c_bound(i));          %me
        pos = inds(r_bound(i)-1,c_bound(i));        %up
    end

    yp_bound(i,1) = pos; 
    yn_bound(i,1) = neg;
    
    if (c_bound(i)-1)<1 || ~face(r_bound(i),c_bound(i)-1)
        %if up is not accesable
        neg = inds(r_bound(i),c_bound(i)+1);          %down
        pos = inds(r_bound(i),c_bound(i));          %me
        % else subtract current entry from next one
    else
        neg = inds(r_bound(i),c_bound(i));          %me
        pos = inds(r_bound(i),c_bound(i)-1);        %up
    end

    xp_bound(i,1) = pos; 
    xn_bound(i,1) = neg;
end





n_in_zs = numel(inface_inds);
sz =3; dev = 2;
gauss = fspecial('gaussian',sz,dev);
% gauss = [-1 -1  -1; -1 8 -1; -1 -1 -1];
% gauss = [0 -1  0; -1 4 0; 0 -1 0];
rhs_reg_mat = lambda1*(z_ref - conv2(z_ref,gauss,'same'));


% rhs_reg_mat = lambda1*(conv2(z_ref,gauss,'same'));
rhs_reg = rhs_reg_mat(inface_inds);
f_w = floor(sz/2);
[boxc, boxr] = meshgrid(-f_w:f_w,-f_w:f_w);

gauss = -gauss;
gauss(2,2) = 1+gauss(2,2);



gaussVec = lambda1*gauss(:);
% gaussVec = -lambda1*gauss(:);
iz_reg = [];
for i=1:n_in_zs
     elems3x3 = inds(sub2ind(size(face),boxr(:)+r_inface(i),boxc(:)+c_inface(i)));
    iz_reg(end+1,:) = elems3x3;
end

nR = 2*sum(in_face(:)) + sum(b_out_full(:));
nC = sum(face(:));
nOnes = sum(in_face(:))*(3+9) + sum(b_out_full(:))*4;
jacobianPattern = sparse([],[],[],nR,nC,nOnes);
constNumber = repmat(1:numel(xp),4,1)';
jacobianPattern(sub2ind([nR nC],constNumber,[xp yp xn yn])) = 1;

offset = numel(yn);
constNumber = repmat(1:numel(xp_bound),4,1)' + offset;
jacobianPattern(sub2ind([nR nC],constNumber,...
    [xp_bound yp_bound xn_bound yn_bound])) = 1;

offset = offset + numel(xp_bound);
constNumber = repmat(1:size(iz_reg,1),9,1)' + offset;
jacobianPattern(sub2ind([nR nC],constNumber,...
    iz_reg)) = 1;



costfun=@(z)depth_cost_nonlin(z,[xp xn],[yp yn],...
[xp_bound xn_bound],[yp_bound yn_bound],...
ncx,ncy,iz_reg,...
im(in_face),rhs_reg,sh_coeff,alb_ref(in_face),gaussVec);

% [cc,pp,qq]=depth_cost_nonlin(z_ref(face),[xp xn],[yp yn],...
%     [xp_bound xn_bound],[yp_bound yn_bound],...
%     ncx,ncy,iz_reg,...
%     im(in_face),rhs_reg,sh_coeff,alb_ref(in_face),gaussVec);
init_z = z_ref(face);
options = optimset('Display','iter-detailed','maxIter',100,'JacobPattern',jacobianPattern);
[z,fval]=lsqnonlin(costfun,init_z,[],[],options);
% zend = z(end);
% z = z(1:end-1)/z(end);
depth = NaN(size(N_ref_in));
depth(face_inds) = z;
offset = mean(depth(inface_inds)) - mean(z_ref(inface_inds));
depth = depth-offset;
depth(~in_face) = NaN;
% % figure; surf(depth,'edgealpha',0);
% 
% if nargin>7
%     z2 = z_gnd(face_inds);
% else
%     z2 = z_ref(face_inds);
% end
% err_ref = abs(A*z2-rhs);
% err_est = abs(A*z-rhs);
% % depth2 = depth;
% % depth2(face_inds) = err_est(1:numel(face_inds));
% % depth22 = depth;
% % depth22(face_inds) = err_ref(1:numel(face_inds));
% s = 1; e = n_zs-n_boundary;
% bnd_s = e+1;
% bnd_e = e+n_boundary-numel(bad_boundary);
% reg_s = bnd_e+1;
% reg_e = bnd_e+n_in_zs;
% if talk
%     figure;plot(s:e+1,err_ref(s:e+1),bnd_s:bnd_e+1,err_ref(bnd_s:bnd_e+1),reg_s:reg_e,err_ref(reg_s:reg_e));title('error in ref depth')
%     fprintf('Ground truth error = %d\n',sum(err_ref.^2));
%     figure;plot(s:e+1,err_est(s:e+1),bnd_s:bnd_e+1,err_est(bnd_s:bnd_e+1),reg_s:reg_e,err_est(reg_s:reg_e));title('error in computed depth')
%     fprintf('Estimated error = %d\n',sum(err_est.^2));
%     face_cost_ref = NaN(size(face));
%     face_cost_ref(constraint_inds) = err_ref(s:e);
%     figure;imagesc(face_cost_ref);title('constraint cost, ref face')
%     face_cost_est = NaN(size(face));
%     face_cost_est (constraint_inds) = err_est(s:e);
%     figure;imagesc(face_cost_est);title('constraint cost, ref est')
% end
% figure;imagesc(depth);
% figure;imagesc(depth2)
