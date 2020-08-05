function active_figure = plot_datadist(data,param,mu_history)
%PLOT_DATADIST Summary of this function goes here
%   Detailed explanation goes here
active_figure =figure();
ax=gca;
axis equal;
hold on;

plot(data(1,:),data(2,:),'b.','color','c');
% set('xlim',[-3 4])
% set('ylim',[0 4.5])
xlim([-3 4]);
ylim([0 4.5]);
xlabel('OFFSET');
ylabel('POWER');

t=-pi:0.01:pi;

for j=1:3
mu=param.mu(j,:);
sigma=param.Sigma(:,(2*j-1):(2*j));
[eigvect,eigval]=eig(sigma);
[tmp,maxpos]=max(mu);

% phi=atan(eigvect(1,maxpos)/eigvect(2,maxpos))
phi=atan(eigvect(2,1)/eigvect(1,1));

x=sqrt(eigval(1,1))*cos(t);
x(2,:)=sqrt(eigval(2,2))*sin(t);
trans=[cos(phi),-sin(phi);sin(phi),cos(phi)];
x=trans*x;
x(1,:)=x(1,:)+mu(1);
x(2,:)=x(2,:)+mu(2);
plot(x(1,:),x(2,:),'LineWidth',2);
ax.ColorOrderIndex=ax.ColorOrderIndex-1;
plot(reshape(mu_history(j,1,:),1,[]),reshape(mu_history(j,2,:),1,[]),'Marker','x','LineWidth',1);
end

hold off;


end

