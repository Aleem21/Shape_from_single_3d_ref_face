function [ pts,tri,rgb,x,y,z,spherical,im_mean,eyes_rgb ] = read_USF_eko( echo_path,r,c,talk )
%READ_USF_EKO path can be an echo file path or a whole directory for batch
if nargin<4
    talk = 0;
end
% f = figure;
if strcmp(echo_path(end-3:end),'.eko')
    [spherical_mean,r,c] = find_USF_spherical(echo_path);
    im_mean = find_USF_im(echo_path);
else
    try
        if echo_path(end)~='/' || echo_path(end)~='\'
            echo_path(end+1)='\';
        end
        echo_path_in = echo_path;
        im_mean = imread([echo_path 'ref_albedo.bmp']);
        load([echo_path_in 'ref_depth.mat']);
    catch
        count = 1;
        echo_paths = getAllFiles(echo_path,1);
        spherical_mean = zeros(r,c);
        im_mean = uint8(zeros(r,c,3));
        for i=1:numel(echo_paths)
            echo_path = echo_paths{i};
            if ~strcmp(echo_path(end-3:end),'.eko')
                continue;
            end
            prev_num = num2str(str2num(echo_path(end-5:end-4))-1);
            prev_name = [echo_path(1:end-6) prev_num '.eko'];
            if (i-4)<1
                continue
            end
            if ~strcmp(echo_paths{i-4},prev_name)
                continue
            end
            [spherical,r,c] = find_USF_spherical(echo_path);
            spherical_mean = spherical_mean + (spherical-spherical_mean)/count;
            im = find_USF_im(echo_path);
            im_mean = im_mean + (im-im_mean)/count;
            count = count + 1;
            %         figure(f);
            %         imshow(im_mean)
            %         pause(0.0001);
        end
        if talk
            fprintf('Models used: %d\n', count-1)
        end
        imwrite(im_mean,[echo_path_in 'ref_albedo.bmp']);
        save([echo_path_in 'ref_depth'],'spherical_mean');
    end
    spherical = [];
    
end

%% convert spherical to cartesian coordinated
[x,y,z] = spherical_to_cart_USF( spherical_mean,r,c );
%% clip to valid region only
valid_region = im2double(imread('D:\Research UCSD\Ravi\Sony SFS\datasets\USF 3D Face Data\USF Raw 3D Face Data Set\data_files\test\range2.bmp'));
valid_region(valid_region==0) = NaN;
im_mean = im2double(im_mean).*repmat(double(valid_region(end:-1:1,:)),1,1,3);

%% plot results
if talk
    figure;surf(x,y,z,im_mean,'FaceColor','interp','edgealpha',0.0);
    axis equal
end

%% output conditioning
offset = min(x(:));
pts = [ y(:)/0.6e4 z(:)/0.6e4 (x(:)-offset)/0.6e4 ]';
x = pts(1,:);
y = pts(2,:);
z = pts(3,:);


%% generate triangulation and rgb
eyes_rgb = (im2double(imread('D:\Research UCSD\Ravi\Sony SFS\datasets\USF 3D Face Data\USF Raw 3D Face Data Set\data_files\test\eyes_mask.bmp')));
eyes_rgb = (reshape(eyes_rgb,size(pts,2),3)');
pts = [x;y;z];
rgb = (reshape(im_mean,size(pts,2),3)');

tri = generate_tri_USF(1:size(im_mean,1),1:size(im_mean,2));
tri = tri(:,:)'-1;

end