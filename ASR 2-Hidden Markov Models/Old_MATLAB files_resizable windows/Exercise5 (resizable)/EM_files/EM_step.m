function [param] = EM_step(data,init_param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% s=init_param.Sigma
% m=init_param.mu
post = posterior(init_param,data');
priors = (ones(1,length(post))*post);

mu=(data(1,:)*post)./priors;
mu=mu';
mu(:,2)=(data(2,:)*post)./priors;

for i = 1:init_param.NumComponents
data_mu{i}=data-mu(i,:)'*ones(1,length(data));

sigmacell{i}=((post(:,i)*ones(1,2))'.*data_mu{i})*data_mu{i}'/priors(i);
end

sigma=sigmacell{1};

for i = 2:init_param.NumComponents
sigma(:,:,i)=sigmacell{i};
end

mixing = priors/size(data,2);

param=gmdistribution(mu,sigma,mixing);
% 
% h = ezcontour(@(x,y)pdf(param,[x y]),[-8 6],[-8 6])
% s=param.Sigma
% m=param.mu

end

