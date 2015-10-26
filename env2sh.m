function [ F_sh_coeff,F_recon ] = env2sh( env,order,talk )
%ENV2SH Summary of this function goes here
%   Detailed explanation goes here
if nargin < 2
    order = 3;
end
if nargin<3
    talk = 0;
end
if size(env,3)>1
    for i=1:size(env,3)
        if nargout >1
            [F_sh_coeff(:,i),F_recon(:,:,i)] = env2sh(env(:,:,i),order,0);
        else
            F_sh_coeff(:,i) = env2sh(env(:,:,i),order,0);
        end
    end
    if talk
        figure;subplot(2,1,1);imshow(uint8(F_recon))
        subplot(2,1,2);imshow(uint8(env))
    end
else
    tdesign = 1;
    regular = 0;
    fliege = 0;
    aziElev2aziIncl = @(dirs) [dirs(:,1) pi/2-dirs(:,2)];
    
    [lon, lat] = meshgrid(linspace(-pi,pi,size(env,2))', linspace(0,pi,size(env,1))');
    sampling_grid = [lon(:) lat(:)];
    
    if tdesign
        [~, tdesign_grid] = getTdesign(2*order);
        tdesign_grid = aziElev2aziIncl(tdesign_grid); % convert to azi-incl
    elseif regular
        regular_grid = grid2dirs(30,30);
        regular_weights = getVoronoiWeights(aziElev2aziIncl(regular_grid));
        
        %         regular_weights = getVoronoiWeights(regular_grid);
    elseif fliege
        [~, fliege_grid, fliege_weights] = getFliegeNodes(order+1);
        fliege_grid = aziElev2aziIncl(fliege_grid); % convert to azi-incl
    end
    
    
    if tdesign
        F = interp2(lon,lat,double(env),tdesign_grid(:,1),tdesign_grid(:,2),'cubic');
        F_sh_coeff = directSHT(order, F, tdesign_grid, 'real', []);
    elseif regular
        F = interp2(lon,lat,double(env),regular_grid(:,1)-pi,regular_grid(:,2),'cubic');
        F_sh_coeff = leastSquaresSHT(5, F, regular_grid, 'real', regular_weights);
    elseif fliege
        F = interp2(lon,lat,double(env),fliege_grid(:,1),fliege_grid(:,2),'cubic');
        F_sh_coeff = directSHT(order, F, fliege_grid, 'real', fliege_weights);
    end
    if nargout >1
        env_sh = inverseSHT(F_sh_coeff, sampling_grid, 'real');
        F_recon = reshape(env_sh,size(env));
        if talk
            figure;subplot(2,1,1);imshow(uint8(F_recon))
            subplot(2,1,2);imshow(uint8(env))
        end
    end
end
end
