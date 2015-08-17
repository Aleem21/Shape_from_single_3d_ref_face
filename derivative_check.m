function [ derivative,max_err,rc,errs] = derivative_check( val,fun,nOnes,type,jackPattern )
%DERIVATIVE_CHECK Summary of this function goes here
%   Detailed explanation goes here
c = numel(val);
cost0 = fun(val);
r = numel(cost0);
thresh = 0;
if nargin<3
    nOnes = r*c;
end
delta = 1e-5;
counter = 1;
rz = zeros(1,nOnes);
cz = zeros(1,nOnes);
vz = zeros(1,nOnes);
for i=1:c
    if mod(i,1000)==0
        fprintf('%d out of %d\n',i,c);
    end
    cur_val = val;
    cur_val(i) = cur_val(i) + delta;
    cost1 = fun(cur_val);
    cur_der = (cost1-cost0)/delta;
    inds = find(abs(cur_der) > thresh);
    vals = cur_der(inds);
    rz(counter:counter+numel(inds)-1) = inds;
    cz(counter:counter+numel(inds)-1) = ones(1,numel(inds))*i;
    vz(counter:counter+numel(inds)-1) = vals;
    counter = counter + numel(inds);
end
rz(counter:end) = [];
cz(counter:end) = [];
vz(counter:end) = [];
derivative = sparse(rz,cz,vz);

if type==2
    [~,d2] = fun(val);
end
max_err = 0;
for i=1:numel(rz)
        if abs(derivative(rz(i),cz(i))-d2(rz(i),cz(i)))>max_err
            max_err = abs(derivative(rz(i),cz(i))-d2(rz(i),cz(i)));
            maxi = rz(i);
            maxj = cz(i);
        end   
end
rc = [maxi; maxj];
errs = [];
if nargin>4
    errs = (jackPattern - (derivative>0))<0;
    errs = sum(errs(:));
end
end

