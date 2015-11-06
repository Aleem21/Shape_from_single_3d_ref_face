% tic
% for i=1:100
%     r = rand(10000,2);
%     for j=1:10
%     getSH(10,r(1:1000,:),'real');
%     end
% end
% toc
tic
for i=1:100
    r = rand(10000,2);
    getSH(10,r,'real');
end
toc
% tic
% for i=1:100
%     costfn_opti(params);
%     
% end
% toc
% t=[];
% for i=1:30
%     dir = rand(i*100,2);
%     tic;
%     for j=1:5
%         getSH(10,dir,'real');
%     end
%     temp = toc;
%     t(end+1) = temp;
% end