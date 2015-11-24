function [ cost,J ] = costfn_skin_model(params, z_ref,chrom_p,D,S,w_o,...
    labels3,wr,wc,alb_mean,alb_var,L_sh_mean,L_sh_var,n_mean,n_var,...
    rho_mean,rho_var,iz_diff,sz,is_face,is_face3,r,type,inds3x3,options )
%COSTFN_SKIN_MODEL Summary of this function goes here
%   Detailed explanation goes here
try
    z_v = params(1:sz(1));
    z2 = nan(size(is_face));
    z2(is_face) = z_v;
    alb_v = params((sz(1)+1):4*sz(1));
    alb = zeros(size(labels3));
    alb(is_face3) = alb_v(:);
    theta = params(4*sz(1)+1:6*sz(1));
    P = params(end);
    L_sh = params(6*sz(1)+1:6*sz(1)+sz(2));
    L_sh = L_sh*P;
    chrom  = params(6*sz(1)+sz(2)+1:6*sz(1)+sz(2)+3);
    delta_a = params(end-6:end-4);
    delta_s1 = params(end-3:end-1);
    n = theta(1:end/2);
    rho_s = theta(end/2+1:end);
    N = normal_from_depth(z2);
    N_v = reshape(N(is_face3),[],3);
    [~,S_r,D_r,alpha,P_dr,Yi,Ft_o,E_lm,d_S_n,l_i,Y_N]=render_DS(N_v',w_o,n,rho_s,alb,r,L_sh,is_face,delta_a,delta_s1,chrom);
    %% Data cost
    S_r = S_r(is_face);
    D_r = D_r(is_face3);
    cost_data_spec= S_r(:)-S(:);
    % cost_data_spec = cost_data_spec;
    cost_data_diffuse = D_r(:)-D(:);
    % cost_data = sum(cost_data_diffuse) + sum(cost_data_spec.^2);
    %% Albedo cost
    % prior
    
    cost_alb_p = zeros(size(is_face));
    for reg=1:10
        alb_r = reshape(alb(labels3==reg),[],3);
        mean_center = alb_r'-repmat(alb_mean{reg},1,size(alb_r,1));
        cost_alb_p(labels3(:,:,1)==reg) = dot(mean_center,(alb_var{reg}\mean_center));
    end
    cost_alb_p = cost_alb_p(is_face);
    % because later on it will be squared and summed by the optimizer itself.
    cost_alb_p = (cost_alb_p).^0.5;
    
    % smoothness
    % right neighbour
    alb1 = reshape(alb(:,1:end-1,:)-alb(:,2:end,:),[],3);
    % neg_cost1 = sum( max(-reshape(alb(:,1:end-1,:),[],3),0),2) + sum( max(-reshape(alb(:,2:end,:),[],3),0),2);
    alb1 = sum(alb1'.^2).*wc{1}.*wr{1}';
    alb1 = [reshape(alb1,size(alb,1),size(alb,2)-1) zeros(size(alb,1),1) ];
    alb1 = alb1(is_face);
    % left neighbour
    % alb2 = reshape(alb(:,2:end,:)-alb(:,1:end-1,:),[],3);
    % alb2 = sum(alb2'.^2).*wc{1}.*wr{1}';
    % alb2 = [reshape(alb2,size(alb,1),size(alb,2)-1) zeros(size(alb,1),1) ];
    % alb2 = alb2(is_face);
    alb2 = [];
    % bot neighbour
    alb3 = reshape(alb(1:end-1,:,:)-alb(2:end,:,:),[],3);
    % neg_cost3 = sum( max(-reshape(alb(1:end-1,:,:),[],3),0),2) + sum( max(-reshape(alb(2:end,:,:),[],3),0),2);
    alb3 = sum(alb3'.^2).*wc{2}.*wr{2}';
    alb3 = [reshape(alb3,size(alb,1)-1,size(alb,2)); zeros(1,size(alb,2))];
    alb3 = alb3(is_face);
    
    % top neighbour
    % alb4 = reshape(alb(2:end,:,:)-alb(1:end-1,:,:),[],3);
    % alb4 = sum(alb4'.^2).*wc{2}.*wr{2}' ;
    % alb4 = [zeros(1,size(alb,2)); reshape(alb4,size(alb,1)-1,size(alb,2))];
    % alb4 = alb4(is_face);
    alb4 = [];
    cost_alb_s = [alb1; alb2; alb3; alb4].^0.5;
    % alb_v_r = alb_v(1:end/3);
    % alb_v_g = alb_v(end/3+1:2*end/3);
    % alb_v_b = alb_v(2*end/3+1:end);
    % cost_alb_s = sum(sum((repmat(alb_v(iz_diff(10,:))',4,1)-alb_v(iz_diff([6 9 11 14],:))).^2));
    
    
    %% Blinn Phong smoothness costs skipped for now
    % smoothness
    cost_rho_s = repmat(rho_s(iz_diff(10,:))',4,1)-rho_s(iz_diff([6 9 11 14],:));
    cost_rho_s = cost_rho_s(:);
    cost_n_s = repmat(n(iz_diff(10,:))',2,1)-n(iz_diff([11 14],:));
    cost_n_s =cost_n_s(:) ;
    % prior
    cost_n_p = ((n-n_mean).^2./n_var).^0.5;
    %% Geometry cost
    % prior
    cost_geo_p = z_v-z_ref;
    % smoothness
    cost_geo_s = repmat(z_v(iz_diff(10,:))',2,1)-z_v(iz_diff([11 14],:));
    cost_geo_s =cost_geo_s(:);
    %% BSSRDF cost
    % absorption parameters
    sc = -[0.0305 0.1236 0.1640]';
    ic = [0.1 0.331 0.605]';
    cost_c = sqrt(sum((delta_a-sc.*delta_s1-ic).^2));
    %% Illumination cost
    cost_illum_p = chrom-chrom_p;
    L_c = L_sh/P-L_sh_mean;
    cost_illum_sh = (L_c'/L_sh_var*L_c).^0.5;
    cost_illum_n = (norm(chrom)-1)*0 + (norm(L_sh/P)-1) + sum(sum(max(-sh2sph(L_sh/P,[50 50]),0)))/5*0 ;
    cost_chrom_n = (norm(chrom)-1);
    
    %% Albedo normalization
    % cost_alb_n = sum((alb(is_face3)).^2)/sz(1)*50 - 1;
    % cost_alb_n = norm(alb(is_face3))/sz(1)*50 - 1;
    cost_alb_n = alb(is_face3)*10/sz(1) - 1;
    
    
    lamb_fp = sqrt(0);
    lamb_hs = sqrt(0);
    lamb_zs = sqrt(0);
    lamb_ln = sqrt(0);
    lamb_fs = sqrt(0);
    lamb_lp = sqrt(0);
    lamb_lsh = sqrt(0);
    lamb_hp = sqrt(0);
    lamb_hc = sqrt(0);
    lamb_zp = sqrt(0);
    
    % lamb_fp = sqrt(10);
    % lamb_fp = sqrt(0.0001);
    % lamb_hs = sqrt(1);
    % lamb_zs = sqrt(0.01);
    % lamb_ln = sqrt(10);
    % lamb_fs = sqrt(10);
    % lamb_lp = sqrt(5);
    % lamb_lsh = sqrt(0.00001);
    % lamb_hp = sqrt(0.5);
    % lamb_hc = sqrt(10);
    % lamb_zp = sqrt(0.0001);
    lamb_fp = sqrt(1);
    lamb_hs = sqrt(100)*0;
    lamb_ln = sqrt(10000);
    lamb_fs = sqrt(1000);
    lamb_lp = sqrt(0.01)*0;
    lamb_lsh = sqrt(1);
    lamb_hp = sqrt(0.5)*0;
    lamb_hc = sqrt(0.1)*0;
    lamb_fn = sqrt(100)*0;
    lamb_diff = 5;
    lamb_spec = 0;
    if type
        lamb_zp = sqrt(0.01);
        lamb_zs = sqrt(0.01);
    else
        lamb_zp = 0;
        lamb_zs = 0;
    end
    cost = [cost_data_diffuse*lamb_diff; cost_data_spec*lamb_spec ;sqrt(sum((cost_alb_p*lamb_fp).^2));...
        cost_alb_s*lamb_fs; cost_n_p*lamb_hp; cost_n_s*lamb_hs;...
        cost_c*lamb_hc; cost_geo_p*lamb_zp; cost_geo_s*lamb_zs; ...
        cost_illum_p*lamb_lp; cost_illum_sh*lamb_lsh; cost_illum_n*lamb_ln;...
        cost_chrom_n*lamb_ln];
%     cost = [cost_data_diffuse*lamb_diff; cost_data_spec*lamb_spec ;cost_alb_p*lamb_fp;...
%         sqrt(sum((cost_alb_s*lamb_fs).^2)); cost_n_p*lamb_hp; cost_n_s*lamb_hs;...
%         cost_c*lamb_hc; cost_geo_p*lamb_zp; cost_geo_s*lamb_zs; ...
%         cost_illum_p*lamb_lp; cost_illum_sh*lamb_lsh; cost_illum_n*lamb_ln;...
%         cost_chrom_n*lamb_ln; cost_alb_n*lamb_fn];
    
    if sum(~isreal(cost))>0
        disp(1)
    end
    
    if nargout >1
        cost_fn_jac = @(params)costfn_skin_model_2(params,chrom_p,D,S,w_o,labels3,...
            sz,is_face,is_face3,r,lamb_hc,lamb_lp,lamb_diff,lamb_spec,type);
        [J_num,ncalls] = compute_Jacobian(params,cost_fn_jac,options);
        % d(diffuse)/d(albedo)
        c_d_alb1 = reshape(inds3x3,[],1);
        P_dr_vec = reshape(repmat(P_dr{1}(:),1,size(inds3x3,2)),[],1);
        P_dg_vec = reshape(repmat(P_dr{2}(:),1,size(inds3x3,2)),[],1);
        P_db_vec = reshape(repmat(P_dr{3}(:),1,size(inds3x3,2)),[],1);
        P_dr_vec(isnan(c_d_alb1)) = [];
        P_dg_vec(isnan(c_d_alb1)) = [];
        P_db_vec(isnan(c_d_alb1)) = [];
        %     v_d_alb = alpha(c_d_alb).*lamb_diff*4;
        %
        %     r_dr_alb = r_d_alb;
        %     r_dg_alb = r_d_alb+sz(1);
        %     r_db_alb = r_d_alb+2*sz(1);
        %     c_dr_alb = c_d_alb+sz(1);
        %     c_dg_alb = c_d_alb+2*sz(1);
        %     c_db_alb = c_d_alb+3*sz(1);
        %     v_dr_alb = v_d_alb*chrom(1).*P_dr_vec;
        %     v_dg_alb = v_d_alb*chrom(2).*P_dg_vec;
        %     v_db_alb = v_d_alb*chrom(3).*P_db_vec;
        r_d_alb = 1:sz(1)*3;
        c_d_alb = sz(1) + (1:sz(1)*3);
        v_d_alb = D_r./alb(is_face3)*lamb_diff;
        % d(diffuse)/d(L_sh_coeff)
        L_sh_offset = sz(1)*6;
        r_dr_L = reshape(repmat(1:sz(1),numel(L_sh),1),[],1);
        c_dr_L = reshape(repmat((1:numel(L_sh))',1,sz(1)),[],1)+L_sh_offset;
        Ft_ov = Ft_o(is_face);
        [row,col] = find(~isnan(inds3x3));
        
        inds = inds3x3(:);
        inds(isnan(inds)) = [];
        r_d_L = repmat(col',numel(L_sh),1);
        c_d_L = repmat((1:numel(L_sh))',1,numel(inds))+sz(1)*6;
        v_dr_L = 4*(Yi(inds,:).* repmat(Ft_ov(inds).*P_dr_vec,1,numel(L_sh)).*repmat(E_lm',numel(inds),1)*chrom(1)*r^2*P*lamb_diff)'.*repmat(alb_v(inds)',numel(L_sh),1);
        v_dg_L = 4*(Yi(inds,:).* repmat(Ft_ov(inds).*P_dg_vec,1,numel(L_sh)).*repmat(E_lm',numel(inds),1)*chrom(2)*r^2*P*lamb_diff)'.*repmat(alb_v(sz(1)+inds)',numel(L_sh),1);
        v_db_L = 4*(Yi(inds,:).* repmat(Ft_ov(inds).*P_db_vec,1,numel(L_sh)).*repmat(E_lm',numel(inds),1)*chrom(3)*r^2*P*lamb_diff)'.*repmat(alb_v(sz(1)*2+inds)',numel(L_sh),1);
        
        % d(diffuse)/d(chrom)
        r_d_chrom = [(1:sz(1))' (1:sz(1))'+sz(1) (1:sz(1))'+sz(1)*2 ];
        c_d_chrom = [ones(sz(1),1) ones(sz(1),1)+1 ones(sz(1),1)+2]+sz(1)*6+numel(L_sh);
        v_d_chrom = lamb_diff*[D_r(1:sz(1))/chrom(1) D_r(sz(1)+1:2*sz(1))/chrom(2) D_r(2*sz(1)+1:3*sz(1))/chrom(3) ];
        v_d_chrom = v_d_chrom*0;
        %d(diffuse)/d(P)
        r_d_P = 1:sz(1)*3;
        c_d_P = ones(sz(1)*3,1)*numel(params);
        v_d_P = lamb_diff*D_r/P*0;
        
        %d(specular)/d(n)
        r_s_n = sz(1)*3+(1:sz(1));
        c_s_n = sz(1)*4+(1:sz(1));
        v_s_n = d_S_n*0;
        
        %d(specular)/d(rho)
        r_s_rho = r_s_n;
        c_s_rho = c_s_n + sz(1);
        v_s_rho = S_r./rho_s*0;
        
        %d(specular)/d(L_sh_coeff)
        r_s_L = repmat((1:sz(1))'+sz(1)*3,1,numel(L_sh));
        c_s_L = repmat((1:numel(L_sh))+sz(1)*6,sz(1),1);
        v_s_L = repmat(rho_s,1,numel(L_sh)).*exp(-l_i'.^2./repmat(n,1,numel(L_sh))/2).*Y_N*P*0;
        
        %d(specular)/d(P)
        r_s_P = (1:sz(1)) + sz(1)*3;
        c_s_P = ones(sz(1),1)*numel(params);
        v_s_P = S_r/P*0;
        
        %d(alb_p)/d(alb)
        v_albp_alb1 = zeros(size(is_face));
        v_albp_alb2 = zeros(size(is_face));
        v_albp_alb3 = zeros(size(is_face));
        for reg=1:10
            alb_r = reshape(alb(labels3==reg),[],3);
            mean_center = alb_r'-repmat(alb_mean{reg},1,size(alb_r,1));
            var_inv = inv(alb_var{reg});
            temp =  mean_center'*[var_inv(1,:)' var_inv(:,1)...
                var_inv(2,:)' var_inv(:,2)...
                var_inv(3,:)' var_inv(:,3)];
            v_albp_alb1(labels3(:,:,1)==reg) = sum(temp(:,1:2),2);
            v_albp_alb2(labels3(:,:,1)==reg) = sum(temp(:,3:4),2);
            v_albp_alb3(labels3(:,:,1)==reg) = sum(temp(:,5:6),2);
        end
        v_albp_alb1 = v_albp_alb1(is_face)./cost_alb_p/2*lamb_fp;
        v_albp_alb2 = v_albp_alb2(is_face)./cost_alb_p/2*lamb_fp;
        v_albp_alb3 = v_albp_alb3(is_face)./cost_alb_p/2*lamb_fp;
        v_albp_alb1(isnan(v_albp_alb1) | isinf(v_albp_alb1) ) = 0;
        v_albp_alb2(isnan(v_albp_alb2) | isinf(v_albp_alb1) ) = 0;
        v_albp_alb3(isnan(v_albp_alb3) | isinf(v_albp_alb1) ) = 0;
        
        r_albp_alb = (1:sz(1))+sz(1)*4;
        c_albp1_alb = (1:sz(1))+sz(1)*1;
        c_albp2_alb = (1:sz(1))+sz(1)*2;
        c_albp3_alb = (1:sz(1))+sz(1)*3;
        
        cost_alb_p_2 = sqrt(sum((cost_alb_p*lamb_fp).^2));
        if cost_alb_p_2>0
            v_albp_alb1_new = (cost_alb_p*lamb_fp .* v_albp_alb1)/cost_alb_p_2;
            v_albp_alb2_new = (cost_alb_p*lamb_fp .* v_albp_alb2)/cost_alb_p_2;
            v_albp_alb3_new = (cost_alb_p*lamb_fp .* v_albp_alb3)/cost_alb_p_2;
            v_albp_alb_new = [v_albp_alb1_new; v_albp_alb2_new; v_albp_alb3_new];
        else
            v_albp_alb_new = zeros(size(cost_alb_p,1)*3,1);
        end
        r_albp_alb_new = r_s_P(end) + ones(size(v_albp_alb_new));
        c_albp_alb_new = sz(1) + (1:sz(1)*3);
        % d(albs)/d(alb)
        diffr = [alb(:,1:end-1,1)-alb(:,2:end,1) zeros(size(is_face,1),1)];
        diffr = diffr(is_face);
        diffg = [alb(:,1:end-1,2)-alb(:,2:end,2) zeros(size(is_face,1),1)];
        diffg = diffg(is_face);
        diffb = [alb(:,1:end-1,3)-alb(:,2:end,3) zeros(size(is_face,1),1)];
        diffb = diffb(is_face);
        r_albs = r_albp_alb_new(end)+(1:sz(1));
        c_albs_pos = (1:sz(1)) + sz(1);
        wcnew = [reshape(wc{1}.*wr{1}',size(is_face,1),size(is_face,2)-1) zeros(size(is_face,1),1)];
        val_posr = diffr.*wcnew(is_face)./alb1.^0.5*lamb_fs;
        val_posg = diffg.*wcnew(is_face)./alb1.^0.5*lamb_fs;
        val_posb = diffb.*wcnew(is_face)./alb1.^0.5*lamb_fs;
        
        face_inds = zeros(size(is_face));
        face_inds(is_face) = 1:sum(is_face(:));
        [row,col] = find(is_face);
        face_inds2 = face_inds;
        face_inds2(end+1,end+1) = 0;
        c_albs_neg = face_inds2(sub2ind(size(face_inds2),row,col+1)) + sz(1);
        
        val_posr(c_albs_neg==0) = 0;
        val_posg(c_albs_neg==0) = 0;
        val_posb(c_albs_neg==0) = 0;
        
        c_albs_neg(c_albs_neg==0) = 1;
        val_posr(isnan(val_posr) | isinf(val_posr) ) = 0;
        val_posg(isnan(val_posg) | isinf(val_posg) ) = 0;
        val_posb(isnan(val_posb) | isinf(val_posb) ) = 0;
        
        diffr = [alb(1:end-1,:,1)-alb(2:end,:,1); zeros(1,size(is_face,2))];
        diffr = diffr(is_face);
        diffg = [alb(1:end-1,:,2)-alb(2:end,:,2); zeros(1,size(is_face,2))];
        diffg = diffg(is_face);
        diffb = [alb(1:end-1,:,3)-alb(2:end,:,3); zeros(1,size(is_face,2))];
        diffb = diffb(is_face);
        r_albs2 = r_albs(end)+(1:sz(1));
        c_albs_pos2 = (1:sz(1)) + sz(1);
        wcnew = [reshape(wc{2}.*wr{2}',size(is_face,1)-1,size(is_face,2)); zeros(1,size(is_face,2))];
        val_posr2 = diffr.*wcnew(is_face)./alb3.^0.5*lamb_fs;
        val_posg2 = diffg.*wcnew(is_face)./alb3.^0.5*lamb_fs;
        val_posb2 = diffb.*wcnew(is_face)./alb3.^0.5*lamb_fs;
        
        c_albs_neg2 = face_inds2(sub2ind(size(face_inds2),row+1,col))+ sz(1);
        val_posr2(c_albs_neg2==0) = 0;
        val_posg2(c_albs_neg2==0) = 0;
        val_posb2(c_albs_neg2==0) = 0;
        c_albs_neg2(c_albs_neg2==0) = 1;
        val_posr2(isnan(val_posr2) | isinf(val_posr2) ) = 0;
        val_posg2(isnan(val_posg2) | isinf(val_posg2) ) = 0;
        val_posb2(isnan(val_posb2) | isinf(val_posb2) ) = 0;
        
        
        % d(n_prior)/d(n)
        r_np_n = r_albs2(end)+(1:sz(1));
        c_np_n = (1:sz(1)) + sz(1)*4;
        v_np_n = 1./(n_var).^0.5*lamb_hp;
        
        % d(n_smoothing)/d(n)
        r_ns_n = repmat((1:sz(1)*2) + r_np_n(end),2,1);
        c_ns_n = [iz_diff([10 11 10 14],:)]+sz(1)*4;
        v_ns_n = [ones(1,sz(1));-1*ones(1,sz(1));ones(1,sz(1));-1*ones(1,sz(1))]*lamb_hs;
        
        % d(L_n)/d(L)
        r_Ln_L = ones(numel(L_sh),1) + r_ns_n(end);
        c_Ln_L = (1:numel(L_sh))+sz(1)*6;
        v_Ln_L = 1/norm(L_sh)*L_sh*lamb_ln;
        
       
        
        % d(z_p)/d(z)
        r_zp_dz = (1:sz(1)) + r_Ln_L(end);
        c_zp_dz = 1:sz(1);
        v_zp_dz = ones(size(c_zp_dz))*lamb_zp;
        
        % d(z_s)/d(z)
        r_zs_dz = repmat((1:sz(1)*2) + r_zp_dz(end),2,1);
        c_zs_dz = [iz_diff([10 11 10 14],:)];
        v_zs_dz = [ones(1,sz(1));-1*ones(1,sz(1));ones(1,sz(1));-1*ones(1,sz(1))]*lamb_zs;
        
        % d(chrom_n)/d(L)
        r_chromn_chrom = ones(3,1) + r_zs_dz(end);
        c_chromn_chrom = (1:3)+sz(1)*6+numel(L_sh);
        v_chromn_chrom = 1/norm(chrom)*chrom*lamb_ln;
        
        
        % d(alb_n)/d(alb)
        r_albn_alb = r_zs_dz(end) + ones(1,sz(1)*3);
        c_albn_alb = sz(1) + (1:sz(1)*3);
        v_albn_alb = 2*alb(is_face3)*lamb_fn*50/sz(1);
        v_albn_alb = lamb_fn*10/sz(1)*ones(size(c_albn_alb));
        
        r = [r_d_alb(:); r_d_L(:);...
            r_d_L(:)+sz(1); r_d_L(:)+sz(1)*2; r_d_chrom(:); r_d_P(:); r_s_n(:);...
            r_s_rho(:); r_s_L(:); r_s_P(:); r_albp_alb_new(:);...
            r_albs(:); r_albs(:); r_albs2(:); r_albs2(:); ...
            r_albs(:); r_albs(:); r_albs2(:); r_albs2(:);  r_albs(:);...
            r_albs(:); r_albs2(:); r_albs2(:); r_np_n(:); r_ns_n(:); r_Ln_L(:);...
            r_chromn_chrom(:); r_zs_dz(:); r_zp_dz(:);];
        %         r = [r_d_alb(:); r_d_L(:);...
        %             r_d_L(:)+sz(1); r_d_L(:)+sz(1)*2; r_d_chrom(:); r_d_P(:); r_s_n(:);...
        %             r_s_rho(:); r_s_L(:); r_s_P(:); r_albp_alb(:); r_albp_alb(:);...
        %             r_albp_alb(:); r_albs(:); r_albs(:); r_albs2(:); r_albs2(:); ...
        %             r_albs(:); r_albs(:); r_albs2(:); r_albs2(:);  r_albs(:);...
        %             r_albs(:); r_albs2(:); r_albs2(:); r_np_n(:); r_ns_n(:); r_Ln_L(:);...
        %             r_chromn_chrom(:); r_zs_dz(:); r_zp_dz(:); r_albn_alb(:)];
        c = [c_d_alb(:); c_d_L(:); c_d_L(:);...
            c_d_L(:); c_d_chrom(:); c_d_P(:); c_s_n(:); c_s_rho(:); c_s_L(:);...
            c_s_P(:); c_albp_alb_new(:); c_albs_pos(:); c_albs_neg(:);...
            c_albs_pos2(:); c_albs_neg2(:); c_albs_pos(:)+sz(1); c_albs_neg(:)+sz(1);...
            c_albs_pos2(:)+sz(1); c_albs_neg2(:)+sz(1); c_albs_pos(:)+sz(1)*2; c_albs_neg(:)+sz(1)*2;...
            c_albs_pos2(:)+sz(1)*2; c_albs_neg2(:)+sz(1)*2; c_np_n(:); c_ns_n(:)...
            ; c_Ln_L(:); c_chromn_chrom(:); c_zs_dz(:); c_zp_dz(:)];
%         c = [c_d_alb(:); c_d_L(:); c_d_L(:);...
%             c_d_L(:); c_d_chrom(:); c_d_P(:); c_s_n(:); c_s_rho(:); c_s_L(:);...
%             c_s_P(:); c_albp1_alb(:); c_albp2_alb(:); c_albp3_alb(:); c_albs_pos(:); c_albs_neg(:);...
%             c_albs_pos2(:); c_albs_neg2(:); c_albs_pos(:)+sz(1); c_albs_neg(:)+sz(1);...
%             c_albs_pos2(:)+sz(1); c_albs_neg2(:)+sz(1); c_albs_pos(:)+sz(1)*2; c_albs_neg(:)+sz(1)*2;...
%             c_albs_pos2(:)+sz(1)*2; c_albs_neg2(:)+sz(1)*2; c_np_n(:); c_ns_n(:)...
%             ; c_Ln_L(:); c_chromn_chrom(:); c_zs_dz(:); c_zp_dz(:);c_albn_alb(:)];
        v = [v_d_alb(:); v_dr_L(:);...
            v_dg_L(:); v_db_L(:); v_d_chrom(:); v_d_P; v_s_n(:); v_s_rho(:);...
            v_s_L(:); v_s_P(:); v_albp_alb_new(:);...
            val_posr(:); -val_posr(:); val_posr2(:); -val_posr2(:); val_posg(:);...
            -val_posg(:); val_posg2(:); -val_posg2(:); val_posb(:); -val_posb(:);...
            val_posb2(:); -val_posb2(:); v_np_n(:) ; v_ns_n(:); v_Ln_L(:);...
            v_chromn_chrom(:); v_zs_dz(:); v_zp_dz(:)];
%         v = [v_d_alb(:); v_dr_L(:);...
%             v_dg_L(:); v_db_L(:); v_d_chrom(:); v_d_P; v_s_n(:); v_s_rho(:);...
%             v_s_L(:); v_s_P(:); v_albp_alb1(:); v_albp_alb2(:); v_albp_alb3(:);...
%             val_posr(:); -val_posr(:); val_posr2(:); -val_posr2(:); val_posg(:);...
%             -val_posg(:); val_posg2(:); -val_posg2(:); val_posb(:); -val_posb(:);...
%             val_posb2(:); -val_posb2(:); v_np_n(:) ; v_ns_n(:); v_Ln_L(:);...
%             v_chromn_chrom(:); v_zs_dz(:); v_zp_dz(:); v_albn_alb(:) ];
        %     tic
        %     J = sparse(r,c,v,numel(cost),numel(params),numel(r));
        %     toc
        %     tic
        J = sparse(r,c,v,numel(cost),numel(params));
        %     toc
        if type
            J(1:sz(1)*4,1:sz(1)) = J_num(1:sz(1)*4,1:sz(1));
            %         J(sz(1)*10+1,sz(1)*6 + numel(L_sh) + 3 + (1:6)) = ...
            J_num(sz(1)*4+1,sz(1)*6 + numel(L_sh) + 3 + (1:6));
            J(sz(1)*13+1 + (1:3),sz(1)*6 + numel(L_sh) + (1:3)) = ...
                J_num(sz(1)*4+1+(1:3),sz(1)*6 + numel(L_sh) + (1:3));
        else
            %         J(sz(1)*10+1,sz(1)*6 + numel(L_sh) + 3 + (1:6)) = ...
            J_num(1,sz(1)*6 + numel(L_sh) + 3 + (1:6));
            J(sz(1)*13+1 + (1:3),sz(1)*6 + numel(L_sh) + (1:3)) = ...
                J_num(1+(1:3),sz(1)*6 + numel(L_sh) + (1:3));
        end
        
        %     J(:,sz(1)*6+numel(L_sh)+(1:3))=0;
    end
catch err
    disp('Something went wrong in cost function calculation')
end
% cost = cost'*cost;
% J = sum(J);
end


