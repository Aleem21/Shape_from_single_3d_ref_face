clear variables
close all

talk = 0;
%% generate ref depth map
sigma = 10 ;
scale = 1;
[dmap_ref, n_ref, N_ref] = generate_ref_depthmap_synth_clean(200,sigma,scale,talk);

%% generate ground truth depth map
sigma = 11;
scale = 1;
[dmap_gnd, n_gnd, N_gnd] = generate_ref_depthmap_synth_clean(200,sigma,scale,talk);

%% generate ref albedo
alb_ref = ~isnan(N_ref);

%% render image
l_ren = [0;0.5;0.25;-1]/2;
talk = 1;
im = render_model_noGL( n_gnd, l_ren, alb_ref,talk);


%% estimate lighting
talk = 0;
is_amb = 1;
non_lin = 0;
l_est = estimate_lighting(n_ref, alb_ref, im,4,is_amb,non_lin);

c4 = render_model_noGL(n_ref,l_est,alb_ref,talk);
for i=1:7
    disp(i)
    depth = estimate_depth(N_ref,alb_ref,im,dmap_ref,l_ren,0.1,'laplac',dmap_gnd,talk);
    % figure;
    % subplot(2,2,1);
    % surf(dmap_gnd,'edgealpha',0.5);title('original');
    % subplot(2,2,2);
    % surf(dmap_ref,'edgealpha',0.5);title('reference');
    % subplot(2,2,3)
    % surf(depth,'edgealpha',0.5);title('reconstructed');
    % subplot(2,2,4)
    % surf(depth-dmap_gnd,'edgealpha',0.5);title('magnified error');
    error = (depth-dmap_gnd).^2;
    n = sum(~isnan(error(:)));

    error(isnan(error)) = 0;
    signal_energy = dmap_ref.^2;
    signal_energy(isnan(signal_energy)) = 0;
    epp = sum(error(:))/n*100;
    epps(i) = epp;
    fprintf('Percentage squared error per pixel = %d %%\n',epp);
    
    [ ~,N_ref2 ] = normal_from_depth( depth );
    %     N_ref = (N_ref+N_ref2)/2;
    N_ref = N_ref2;
    depths{i} = depth;
    errs{i} = N_ref-N_gnd;
    error = (N_ref-N_gnd).^2;
    error(isnan(error)) = 0;
    nreferr(i) = sum(error(:))/n*100;
    fprintf('error in Nref = %d\n', sum(error(:))/n*100);
end
figure;
subplot(2,2,1);
surf(dmap_gnd,'edgealpha',0.5);title('original');
subplot(2,2,2);
surf(dmap_ref,'edgealpha',0.5);title('reference');
subplot(2,2,3)
surf(depth,'edgealpha',0.5);title('reconstructed');
subplot(2,2,4)
surf(depth-dmap_gnd,'edgealpha',0.5);title('magnified error');
error = (depth-dmap_gnd).^2;

error(isnan(error)) = 0;
signal_energy = dmap_ref.^2;
signal_energy(isnan(signal_energy)) = 0;
epp = sum(error(:))/sum(signal_energy(:))*100;
fprintf('Percentage squared error per pixel = %d %%\n',epp);

figure;
while(1)
    for i=1:numel(depths)
        surf(errs{i})
        title(num2str(i))
        drawnow
        pause(0.5)
    end
end
