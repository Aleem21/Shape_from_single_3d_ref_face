function [ costfun, face,nData,nBound,nReg,jacobianPattern ] = get_depth_costfun( z_ref, im,alb_ref, sh_coeff, eye_mask,lambda1,lambda_bound,type,is_dz_depth,is_l_sh,lambda_dz,lambda_reg2)
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


if type==2
    ncx = -lambda_bound*ncx;
    ncy = -lambda_bound*ncy;
    ncx(b_in_full) = -ncx_in(b_in_full);
    ncy(b_in_full) = -ncy_in(b_in_full);
else
    ncx = -ncx(b_out_full);
    ncy = -ncy(b_out_full);
end

rho_ref = alb_ref(in_face);

% r,c to index number in face Z's
% [fr,fc] = meshgrid(1:size(face,1),1:size(face,2));
% sub2ind_face = find_elements([r_face c_face],fr',fc');

sub2ind_face = zeros(size(face));
sub2ind_face(face) = 1:sum(face(:));

%% cost function
% Data term
xp = sub2ind_face(sub2ind(size(im),r_inface,max(c_inface-1,1)));
xn = sub2ind_face(sub2ind(size(im),r_inface,c_inface));

yp = sub2ind_face(sub2ind(size(im),max(r_inface-1,1),c_inface));
yn = sub2ind_face(sub2ind(size(im),r_inface,c_inface));

ind = find(xp==0);
for i = 1:numel(ind)
    xp(ind(i))=  xn(ind(i));
    xn(ind(i))=sub2ind_face(sub2ind(size(im),r_face(ind(i)),c_face(ind(i))+1));
end
ind = find(yp==0);
for i = 1:numel(ind)
    yp(ind(i))=  yn(ind(i));
    yn(ind(i))= sub2ind_face(sub2ind(size(im),r_face(ind(i))+1,c_face(ind(i))));
end

% boundary term
[r_bound,c_bound] = find(b_out_full);
num_boundaries = sum(b_out_full(:));
yp_bound = zeros(num_boundaries,1);
yn_bound = zeros(num_boundaries,1);

xp_bound = zeros(num_boundaries,1);
xn_bound = zeros(num_boundaries,1);
for i=1:num_boundaries
    if type==1
        [xp_bound(i), xn_bound(i), yp_bound(i), yn_bound(i) ] ...
            = find_grad_inds_boundary(r_bound(i),c_bound(i),face,sub2ind_face,[-1,-1]);
    else
        r = r_bound(i);
        c = c_bound(i);
        
        [xp_x_y, xn_x_y, yp_x_y, yn_x_y,pref_v ] ...
            = find_grad_inds_boundary(r,c,face,sub2ind_face,[1,1]);
        [xp_xp_y, xn_xp_y, yp_xp_y, yn_xp_y ] ...
            = find_grad_inds_boundary(r,c-pref_v(1),face,sub2ind_face,pref_v);
        
        [xp_x_yp, xn_x_yp, yp_x_yp, yn_x_yp ] ...
            = find_grad_inds_boundary(r-pref_v(2),c,face,sub2ind_face,pref_v);
        nc = [ncx(r,c); ncy(r,c)];
        nc_xp_y = [ncx(r,c-pref_v(1)); ncy(r,c-pref_v(1))]*0+nc;
        nc_x_yp = [ncx(r-pref_v(2),c); ncy(r-pref_v(2),c)]*0+nc;
        %         xn = c-pref_v(2);
        %         yn = r-pref_v(1);
        i_bound(i,:) = [   xp_xp_y     xn_xp_y ...
            yp_xp_y     yn_xp_y ...
            xp_x_y      xn_x_y  ...
            yp_x_y      yn_x_y  ...
            xp_x_yp     xn_x_yp ...
            yp_x_yp     yn_x_yp ...
            xp_x_y      xn_x_y  ...
            yp_x_y      yn_x_y  ...
            ];
        val_bound(i,1:8) = [nc_xp_y(1)*nc(1) -nc_xp_y(1)*nc(1)...
            nc_xp_y(2)*nc(1) -nc_xp_y(2)*nc(1)...
            -nc(1)*nc(1)       nc(1)*nc(1)...
            -nc(2)*nc(1)       nc(2)*nc(1)]*pref_v(1);
        val_bound(i,9:16) =[nc_x_yp(1)*nc(2) -nc_x_yp(1)*nc(2)...
            nc_x_yp(2)*nc(2) -nc_x_yp(2)*nc(2)...
            -nc(1)*nc(2)       nc(1)*nc(2)...
            -nc(2)*nc(2)       nc(2)*nc(2)]*pref_v(2);
        if sum(isnan(val_bound(i,:)))>0
            val_bound(i,:) = 0;
        end
        % if pref_v== 1  -> bellow are the +ve terms
        % if pref_v==-1  -> bellow are the -ve terms
        
        
    end
end


% regularization term
% in_inface = (in_face-b_in_full);
in_inface = in_face;

[r_innface,c_innface] = find(in_inface);
% innface_inds = sub2ind(size(face),r_innface,c_innface);
innface_inds = find(in_inface);

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

if nargout >5
    % Jacobian Pattern
    nR = sum(in_inface(:))+sum(in_inface(:)) + sum(b_out_full(:)) + sum(in_inface(:));
    nC = sum(face(:));
    nOnes = sum(in_face(:))*(3+9+9) + sum(b_out_full(:))*4;
    jacobianPattern = sparse([],[],[],nR,nC,nOnes);
    constNumber = repmat(1:numel(xp),4,1)';
    jacobianPattern(sub2ind([nR nC],constNumber,[xp yp xn yn])) = 1;
    
    offset = numel(yn);
    if type==1
        constNumber = repmat(1:numel(xp_bound),4,1)' + offset;
        jacobianPattern(sub2ind([nR nC],constNumber,...
            [xp_bound yp_bound xn_bound yn_bound])) = 1;
    else
        constNumber = repmat(1:size(val_bound,1),16,1)' + offset;
        jacobianPattern(sub2ind([nR nC],constNumber,...
            i_bound)) = 1;
    end
    
    offset = offset + numel(xp_bound);
    constNumber = repmat(1:size(iz_reg,1),9,1)' + offset;
    jacobianPattern(sub2ind([nR nC],constNumber,...
        iz_reg)) = 1;
    
    
    offset = offset + size(iz_reg,1);
    constNumber = repmat(1:size(iz_reg,1),9,1)' + offset;
    jacobianPattern(sub2ind([nR nC],constNumber,...
        iz_reg)) = 1;
end
% get eye map
if type==1
    costfun=@(z)cost_nonlin_depth(z,[xp xn],[yp yn],...
        [xp_bound xn_bound],[yp_bound yn_bound],...
        ncx,ncy,iz_reg,...
        im(in_face),rhs_reg,lambda_dz,lambda_reg2,sh_coeff,rho_ref,z_ref(in_face),gaussVec,type,is_dz_depth,eye_mask(in_face));
else
    costfun=@(z)cost_nonlin_depth(z,[xp xn],[yp yn],...
    [],[],...
    ncx,ncy,iz_reg,...
    im(in_face),rhs_reg,lambda_dz,lambda_reg2,sh_coeff,rho_ref,z_ref(in_face)...
    ,gaussVec,type,is_dz_depth,is_l_sh,eye_mask(in_face),i_bound,val_bound,in_face,face);

end
nData = numel(xp);
nBound = numel(xp_bound);
nReg = size(iz_reg,1);
end

