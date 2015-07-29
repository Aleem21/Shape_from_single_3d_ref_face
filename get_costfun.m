function [ costfun, face,nData,nBound,nReg,jacobianPattern ] = get_costfun( z_ref, im,alb_ref, sh_coeff, lambda1)
%GET_COSTFUN Summary of this function goes here
%   Detailed explanation goes here

%% Pre processing
[ face,face_inds, inface_inds,in_face,r_face,c_face,r_inface,c_inface,...
    b_out_full,b_in_full ] = preprocess_estimate_depth( z_ref );


% boundary contour normals
[ncy,ncx] = find_countour_normal(b_out_full);
ncx(~b_out_full) = NaN;
ncy(~b_out_full) = NaN;

[ncy_in,ncx_in] = find_countour_normal(b_in_full);
ncx_in(~b_in_full) = NaN;
ncy_in(~b_in_full) = NaN;

ncx = -ncx(b_out_full);
ncy = -ncy(b_out_full);
rho_ref = alb_ref(in_face);

% r,c to index number in face Z's
[fr,fc] = meshgrid(1:size(face,1),1:size(face,2));
sub2ind_face = find_elements([r_face c_face],fr',fc');

%% cost function
% Data term
xp = sub2ind_face(sub2ind(size(im),r_inface,c_inface-1));
xn = sub2ind_face(sub2ind(size(im),r_inface,c_inface));

yp = sub2ind_face(sub2ind(size(im),r_inface-1,c_inface));
yn = sub2ind_face(sub2ind(size(im),r_inface,c_inface));
[r_bound,c_bound] = find(b_out_full);

% boundary term
num_boundaries = sum(b_out_full(:));
yp_bound = zeros(num_boundaries,1);
yn_bound = zeros(num_boundaries,1);

xp_bound = zeros(num_boundaries,1);
xn_bound = zeros(num_boundaries,1);
for i=1:num_boundaries
    if (r_bound(i)-1)<1 || ~face(r_bound(i)-1,c_bound(i))
        %if up is not accesable
        neg = sub2ind_face(r_bound(i)+1,c_bound(i));          %down
        pos = sub2ind_face(r_bound(i),c_bound(i));          %me
        % else subtract current entry from next one
    else
        neg = sub2ind_face(r_bound(i),c_bound(i));          %me
        pos = sub2ind_face(r_bound(i)-1,c_bound(i));        %up
    end
    yp_bound(i,1) = pos;
    yn_bound(i,1) = neg;
    
    if (c_bound(i)-1)<1 || ~face(r_bound(i),c_bound(i)-1)
        %if up is not accesable
        neg = sub2ind_face(r_bound(i),c_bound(i)+1);          %right
        pos = sub2ind_face(r_bound(i),c_bound(i));          %me
        % else subtract current entry from next one
    else
        neg = sub2ind_face(r_bound(i),c_bound(i));          %left
        pos = sub2ind_face(r_bound(i),c_bound(i)-1);        %up
    end
    
    xp_bound(i,1) = pos;
    xn_bound(i,1) = neg;
end


% regularization term
in_inface = (in_face);
[r_innface,c_innface] = find(in_inface);
innface_inds = sub2ind(size(face),r_innface,c_innface);


sz = 3; dev = 2;
gauss = fspecial('gaussian',sz,dev);
rhs_reg_mat = lambda1*(z_ref - conv2(z_ref,gauss,'same'));

rhs_reg = rhs_reg_mat(innface_inds);
f_w = floor(sz/2);
[boxc, boxr] = meshgrid(-f_w:f_w,-f_w:f_w);
modified_gauss = diag([0 1 0])-gauss; % 1- gauss



gaussVec = lambda1*modified_gauss(:);
iz_reg = zeros(numel(r_innface),numel(boxc));
for i=1:numel(r_innface)
    elems3x3 = sub2ind_face(sub2ind(size(face),boxr(:)+r_innface(i),boxc(:)+c_innface(i)));
    iz_reg(i,:) = elems3x3;
end

if nargout >2
    % Jacobian Pattern
    nR = sum(in_face(:))+sum(in_inface(:)) + sum(b_out_full(:));
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
end

costfun=@(z)depth_cost_nonlin(z,[xp xn],[yp yn],...
    [xp_bound xn_bound],[yp_bound yn_bound],...
    ncx,ncy,iz_reg,...
    im(in_face),rhs_reg,sh_coeff,rho_ref,gaussVec);

nData = numel(xp);
nBound = numel(xp_bound);
nReg = size(iz_reg,1);
end

