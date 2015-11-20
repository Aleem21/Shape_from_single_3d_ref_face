function [ cost ] = costfn_skin_model_2(params,chrom_p,D,S,w_o,labels3,...
    sz,is_face,is_face3,r,lamb_hc,lamb_lp,lamb_diff,lamb_spec,type)
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
    [~,S_r,D_r]=render_DS(N_v',w_o,n,rho_s,alb,r,L_sh,is_face,delta_a,delta_s1,chrom);
    %% Data cost
    S_r = S_r(is_face);
    D_r = D_r(is_face3);
    cost_data_spec= S_r(:)-S(:);
    cost_data_diffuse = D_r(:)-D(:);
    
    %% BSSRDF cost
    % absorption parameters
    sc = -[0.0305 0.1236 0.1640]';
    ic = [0.1 0.331 0.605]';
    cost_c = sqrt(sum((delta_a-sc.*delta_s1-ic).^2));
    %% Illumination cost
    cost_illum_p = chrom-chrom_p;
    if type
        cost = [cost_data_diffuse*lamb_diff; cost_data_spec*lamb_spec; cost_c*lamb_hc; ...
            cost_illum_p*lamb_lp;];
    else
        cost = [cost_c*lamb_hc; cost_illum_p*lamb_lp;];
        
    end
    if sum(~isreal(cost))>0
        disp(1)
    end
    
catch err
    disp('Something went wrong in cost function calculation')
end
end


