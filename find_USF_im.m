function [ im ] = find_USF_im( echo_path )
%FIND_USF_IM reads the corresponding image of the eko file in USF database

impath = [echo_path(1:end-3) 'jpg'];
im = imread(impath);
im = im(end:-1:1,:,:);
end

