function [ out,r,c ] = find_USF_spherical( echo_path,talk )
%FIND_USF_SPHERICAL reads a single eko file from USF dataset and gives its
%spherical coordinate radius values


if nargin<2
    talk = 0;
end
fid = fopen(echo_path);
infeed = fread(fid);
fclose(fid);
intext = char(infeed');

%% parse the header
offset = strfind(intext,'DATA=')+5; offset = offset(1);
header = intext(1:offset);
if talk
    fprintf('%s\n',header);
end
r_i = strfind(header,'NLT=')+4;
r_i_e = strfind(header(r_i:end),char(10))+r_i-2;
r = str2num(header(r_i:r_i_e(1)));

c_i = strfind(header,'NLG=')+4;
c_i_e = strfind(header(c_i:end),char(10))+c_i-2;
c = str2num(header(c_i:c_i_e(1)));

elems = r*c;

little = reshape(infeed((2:2:elems*2+1)+offset),r,c);
big = reshape(infeed((1:2:elems*2)+offset),r,c);
big(big==128) = nan;
out = (big*2^8+little);

end

