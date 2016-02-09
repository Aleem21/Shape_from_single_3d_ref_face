function [ im ] = render_model_noGL( n, sh_coeff, alb,talk )
%RENDER_MODEL 
% [ im ] = render_model_noGL( n, sh_coeff, alb,talk )
if nargin<4
    talk = 0;
end
if size(sh_coeff,1)>1
    sh_coeff = sh_coeff';
end
if size(alb,3)>1
    for i=1:size(alb,3)
        im(:,:,i) = render_model_noGL(n,sh_coeff,alb(:,:,i),0);
    end
else
    sh_coeff(10) = 0;
    sh_coeff(10) = [];
    alb = alb(:)';
    nx = n(:,:,1);
    nx = nx(:)';        %change to row vector
    ny = n(:,:,2);
    ny = ny(:)';
    nz = n(:,:,3);
    nz = nz(:)';
    
    Y = [ones(size(nx));
        nx;
        ny;
        nz;
        nx.*ny;
        nx.*nz;
        ny.*nz;
        (nx.^2 - ny.^2);
        (3*nz.^2-1) ];
    sh_coeff1 = sh_coeff(1);
    sh_coeff2 = sh_coeff(2:end);
    Y1 = Y(1,:);
    Y2 = Y(2:end,:);
    im = alb.*(sh_coeff1*Y1 +  sh_coeff2*Y2);
%     im(im<0) = 0;
    im = reshape(im,size(n,1),size(n,2));
end
if talk
    figure; imshow(im);
    title('rendered, NO GL')
end
end

