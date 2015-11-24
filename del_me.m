% optimum = lsqnonlin(costfn_opti,params,[],[],options);

% tic
% J = sparse(1:100000,1,1,100000,100000);
% toc

J = spalloc(10000,10000,100000);
t = zeros(1,5);
tic
tic
for i=1:20
    J((i-1)*5000 + (1:5000),1) = 1;
    t(i)=toc;
end
figure;plot(t)
for i=1:100
    tic
    J = sparse(1:50000,1,1,100000+(i-1)*10000,100000+(i-1)*10000);
    t2(i) = toc;
end
figure;plot(t2)
% [ z_o,alb_o,n_o,rho_o,L_sh_o,chrom_o,delta_a_o,delta_s1_o,P_o ] = ...
%     skin_model_optimization( is_face,L_sh_mean,L_sh_var,z_ref,chrom_p,...
%     alb_init,n_init,rho_init,L_sh_init,D,S,w_o,labels,im,alb_mean_cell,...
%     alb_var,n_mean,n_var,rho_mean,rho_var,delta_a_init,delta_s1_init,P_init,rad,0 );

% running = 1;
% J_small = J(1:numel([cost_data_diffuse*20; cost_data_spec]),:);
% calls = 0;
% 
% % while(calls<3)
% %     calls = calls+1;
% %     disp(calls)
% %     J2 = J_small>0;
% %     vars = zeros(1,size(J_small,1));
% %     rows = zeros(1,size(J_small,1));
% %     tic
% %     ind = 0;
% %     while(nnz(J_small)>0)
% %         ind = ind+1;
% %         [c,r] = find(J2',1,'first');
% %         if isempty(c)
% %             break
% %         end
% %         vars(ind) = c;
% %         rows(ind) = r;
% %
% %         J2(:,J2(r,:)) = 0;
% %     end
% %
% %     vars(ind+1:end)=[];
% %     rows(ind+1:end)=[];
% %     J_small(sub2ind(size(J_small),rows,vars)) = 0;
% %     toc
% %
% J3 = J;
% iter = 1;
% inps{iter} = [];
% outs{iter} = [];
% for i=1:size(J,2)
%     inps{i} = [];
%     outs{i} = [];
%     J2 = J3;
%     for j=i:size(J,2)
%         t =find(J2(:,j));
%         inps{i}(end+1:end+numel(t)) = ones(1,numel(t));
%         outs{i}(end+1:end+numel(t)) = t;
%         dependent = sum(J2(J2(:,j)>0,:))>0;
%         if nnz(dependent)>0
%             disp('')
%         end
%         J2(:,dependent) = 0;
%         if nnz(J2)==0
%             break
%         end
%     end
%     J3(sub2ind(size(J3),outs{i},inps{i})) = 0;
%     if nnz(J3)==0
%         break
%     end
%         
% end
% % end