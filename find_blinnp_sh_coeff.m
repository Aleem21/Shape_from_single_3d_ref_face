function [ sh_coeff ] = find_blinnp_sh_coeff( N,w_o,w_i_grid,n,Nord )
%FIND_BLINNP_SH_COEFF Computes the SH_coefficients for the lumped terms of
%brdf and A ^\prime for a point on surface with normal N.

A = max(cos(w_i_grid(:,2)),0).*(cos(find_h_local(N,w_o,w_i_grid)).^n);
sh_coeff = directSHT(Nord, A, w_i_grid, 'real', []);


end

