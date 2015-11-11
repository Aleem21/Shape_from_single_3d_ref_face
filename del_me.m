running = 1;
J_small = J(1:numel([cost_data_diffuse*20; cost_data_spec]),:);
calls = 0;

% while(calls<3)
%     calls = calls+1;
%     disp(calls)
%     J2 = J_small>0;
%     vars = zeros(1,size(J_small,1));
%     rows = zeros(1,size(J_small,1));
%     tic
%     ind = 0;
%     while(nnz(J_small)>0)
%         ind = ind+1;
%         [c,r] = find(J2',1,'first');
%         if isempty(c)
%             break
%         end
%         vars(ind) = c;
%         rows(ind) = r;
%
%         J2(:,J2(r,:)) = 0;
%     end
%
%     vars(ind+1:end)=[];
%     rows(ind+1:end)=[];
%     J_small(sub2ind(size(J_small),rows,vars)) = 0;
%     toc
%
J3 = J;
iter = 1;
inps{iter} = [];
outs{iter} = [];
for i=1:size(J,2)
    inps{i} = [];
    outs{i} = [];
    J2 = J3;
    for j=i:size(J,2)
        t =find(J2(:,j));
        inps{i}(end+1:end+numel(t)) = ones(1,numel(t));
        outs{i}(end+1:end+numel(t)) = t;
        dependent = sum(J2(J2(:,j)>0,:))>0;
        if nnz(dependent)>0
            disp('')
        end
        J2(:,dependent) = 0;
        if nnz(J2)==0
            break
        end
    end
    J3(sub2ind(size(J3),outs{i},inps{i})) = 0;
    if nnz(J3)==0
        break
    end
        
end
% end