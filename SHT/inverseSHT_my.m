function F = inverseSHT_my(F_N, dirs, basisType,Y_N)
%INVERSE_SHT Perform the inverse spherical harmonic transform
%
%   N:  maximum order of harmonics
%   F: the spherical function recreated at directions 'dirs'
%   dirs:   [azimuth inclination] angles in rads for each evaluation point,
%           where inclination is the polar angle from zenith
%           theta = pi/2-elevation
%   basisType:  'complex' or 'real' spherical harmonics
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 10/10/2013
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin <4
        % vector of spherical harmonics
        N = sqrt(size(F_N,1)) - 1;
        Y_N = getSH(N, dirs, basisType);
    end
    % perform the inverse transform up to degree N
    if size(F_N,2)==1
        F = Y_N*F_N;
    else
        F = dot(Y_N',F_N);
    end

end
