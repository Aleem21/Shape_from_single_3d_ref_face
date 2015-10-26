sz = [50 50];
[lon, lat] = meshgrid(linspace(-pi,pi,sz(2))', linspace(0,pi,sz(1))');
    sampling_grid = [lon(:) lat(:)];
    sampling_grid_sp = aziElev2aziIncl(sampling_grid);
    [N_test(1,:),N_test(2,:),N_test(3,:)] = sph2cart(sampling_grid_sp(:,1),sampling_grid_sp(:,2),ones(size(sampling_grid_sp,1),1));

    
    N_test = N_v;

for i=1:size(N_test,2)
    
    if mod(i,100)==0
        disp(i/size(N_test,2)*100);
    end
    if sum(isnan(N_test(:,i)))
        continue
    end
    N_dot_H = dot(repmat(N_test(:,i),1,size(H_all,2)),H_all);
    A = rho_s*(n+2)/(2*pi)* Light_wi'.*(N_dot_H.^n);
    A_sh = directSHT(Nord, A', tdesign_grid, 'real', []);
    
    
    
    A = max(dot(repmat(N_v(:,i),1,size(wi_all,2)),wi_all),0);
    A_lm = directSHT(Nord, A', tdesign_grid, 'real', []);
%     
%     l = zeros((Nord+1)^2,1);
%     for j=0:Nord
%         l(j^2+1:(j+1)^2) = j;
%     end
%     const = sqrt(4*pi./(2*l+1));
%     
%     E_lm = const.*A_lm;
%     for j=0:Nord
%         E_lm(i^2+1:(j+1)^2) = E_lm(median([j^2+1 (j+1)^2]));
%     end
%     
    
    
    [temp(1) temp(2)]=cart2sph(N_test(1,i),N_test(2,i),N_test(3,i));

    temp = aziElev2aziIncl(temp);
    S(i) = inverseSHT(A_sh.*E_lm,temp,'real');
end
S(size(N_test,2))=0;

S = reshape(S,sz);




Light_wi = inverseSHT(Light_sh, sampling_grid, 'real');
[wix_all,wiy_all,wiz_all] = sph2cart(sampling_grid_sp(:,1),sampling_grid_sp(:,2),1);
wi_all = [wix_all';wiy_all';wiz_all'];
H_all = (repmat(w_o ,1,size(wi_all,2))+wi_all);
H_all = H_all./repmat(sqrt(sum(H_all.^2)),3,1);
i = 100;
[Nix,Niy,Niz] = sph2cart(grid_sph(i,1),grid_sph(i,2),1);
    Ni = [Nix;Niy;Niz];
    N_dot_H = dot(repmat(Ni,1,size(H_all,2)),H_all);
%     N_dot_wi = dot(repmat(Ni,1,size(wi_all,2)),wi_all);
    A = rho_s*(n+2)/(2*pi)* Light_wi'.*(N_dot_H.^n);
%     A = Light_wi';
    A_sh = directSHT(Nord, A', tdesign_grid, 'real', []);



S = inverseSHT(A_sh, rec_pts, 'real');
S = reshape(S,size(N,1),size(N,2));

AA = inverseSHT(A_sh, sampling_grid, 'real');
AA = reshape(AA,sz);
