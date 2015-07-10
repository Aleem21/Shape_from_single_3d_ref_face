clear variables
close all
import plot3D_helper.label_axis

talk = 4;
%% generate ref depth map
type = 'face';
switch type
    case 'face'
        ply_path = '.\data\ref_model.ply';
        cRes = 168*2;
        rRes = 192*2;
        % adjust face model
        Rneutral = makehgtform('translate',[0 0.2 0])* makehgtform('zrotate',deg2rad(-7))*makehgtform('xrotate',deg2rad(20))*makehgtform('yrotate',deg2rad(10));
        RScale = makehgtform('scale',[cRes rRes min([cRes rRes])]/2)*makehgtform('translate',[1 1 -0.5]);
        Rpose = RScale * Rneutral;

        xrange = [0 cRes];
        yrange = [0 rRes];
        [dmap_ref, n_ref, N_ref] = generate_ref_depthmap(ply_path, rRes, cRes, Rpose,[],xrange,yrange, talk);

    case 'sphere'
        ply_path = '.\data\sphere.ply';
        cRes = 20;
        rRes = 20;
        xrange = [-150 150];
        yrange = [-150 150];
        Rpose = eye(4);
        [dmap_ref, n_ref, N_ref] = generate_ref_depthmap(ply_path, rRes, cRes, Rpose,[],xrange,yrange, talk);
    case 'synth'
        [dmap_ref, n_ref, N_ref] = generate_ref_depthmap_synth(100,talk);

    otherwise
        error('no matching case found');
end


%% generate ref albedo
alb_ref = ~isnan(N_ref);

%% render image
l_ren = [0;0.5;0.25;-1]/2;
talk = 1;
GL = false;
type = 'face';
if GL
    switch type
        case 'face'
            ply_path = '.\data\ref_model.ply';
            cRes = 168*2;
            rRes = 192*2;
            % adjust face model
            Rpose = makehgtform('translate',[0 0.2 0])* makehgtform('zrotate',deg2rad(-7))*makehgtform('xrotate',deg2rad(20))*makehgtform('yrotate',deg2rad(10));
            
            xrange = [-1 1];
            yrange = [-1 1];
            im = render_model_general(ply_path, l_ren, Rpose, rRes, cRes, xrange, yrange, talk);
        case 'sphere'
            ply_path = '.\data\sphere.ply';
            cRes = 20;
            rRes = 20;
            xrange = [-150 150];
            yrange = [-150 150];
            Rpose = eye(4);
            im = render_model_general(ply_path, l_ren, Rpose, 200, 200, xrange, yrange, talk);
        otherwise
            error('no matching case found');
    end
else
    im = render_model_noGL( n_ref, l_ren, alb_ref,talk);
end


%% estimate lighting
talk = 0;
l_est = estimate_lighting(n_ref, alb_ref, im,4);
% c4 = render_model('./data/sphere.ply',l_est,talk,Rpose, 150,150);
epp_arr = [];
c4 = render_model_noGL(n_ref,l_est,alb_ref,talk);
std_list = 5;
for std = std_list
    fprintf('Std : %d\n',std);
    noise = randn(size(dmap_ref))*std;
    depth = estimate_depth(N_ref,alb_ref,im,dmap_ref+noise,l_est,1,1);
    
    figure('name',['sigma' num2str(std)]);
    subplot(2,2,1);
    surf(dmap_ref,'edgealpha',0.5);title('original');axis equal
    subplot(2,2,2);
    surf(dmap_ref+noise,'edgealpha',0.5);title('reference');axis equal
    subplot(2,2,3)
    surf(depth,'edgealpha',0.5);title('reconstructed');axis equal
    subplot(2,2,4)
    surf(depth-dmap_ref,'edgealpha',0.5);title('magnified error');axis equal
    error = (depth-dmap_ref).^2;
    
    error(isnan(error)) = 0;
    signal_energy = dmap_ref.^2;
    signal_energy(isnan(signal_energy)) = 0;
    epp = sum(error(:))/sum(signal_energy(:))*100;
    fprintf('Percentage squared error per pixel = %d %%\n',epp);
    fprintf('Singal energy = %d \n',sum(signal_energy(:)));
    fprintf('Noise energy = %d \n',sum(sum(noise.^2)));
    epp_arr(end+1) = epp;
end
figure; plot(std_list,epp_arr)