function image = sh2sph(sh_coeff,sz)


if size(sh_coeff,2)>1
    for i=1:size(sh_coeff,2)
        image(:,:,i) = sh2sph(sh_coeff(:,i),sz);
    end
else
    [lon, lat] = meshgrid(linspace(-pi,pi,sz(2))', linspace(0,pi,sz(1))');
    sampling_grid = [lon(:) lat(:)];
    
    env_sh = inverseSHT(sh_coeff, sampling_grid, 'real');
    image = reshape(env_sh,sz);
end

end