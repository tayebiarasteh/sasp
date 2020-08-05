close all

%init pdfs
mu=[1,1;2,2;3,3];
sigma1=[1,0;0,1];
sigma2=[1,0;0,1];
sigma3=[1,0;0,1];
sigma=sigma1;
sigma(:,:,2)=sigma2;
sigma(:,:,3)=sigma3;

%create common pdf object
param1=gmdistribution(mu,sigma);
param2=gmdistribution(mu,sigma);

%create datapoints
train1= simulate_word1();
train2= simulate_word2();

%array to store mu in every EM-step
% mu_history=[];



for i=1:50
    param1=EM_step(train1,param1);
    param2=EM_step(train2,param2);
end

%classifie training data
[tmp1,tmp2,log1]=bayes_classifier(param1,param2,[0.5,0.5],train1);
[tmp1,tmp2,log2]=bayes_classifier(param1,param2,[0.5,0.5],train2);

x=[-2,4];
y=[0,4];
z=[-2,1];
precission = 0.05;

%plot decission surface
a = calculate_decission(x,y,precission,param1,param2,[0.5,0.5]);
surf(x(1):precission:x(2),y(1):precission:y(2),a,'LineWidth',0.01,'EdgeAlpha',0.4)

hold on
% a = a*40/max(max(abs(a)))+40;
colormap_factor = max(max(abs(a)));

caxis([-colormap_factor, colormap_factor]);

colormap(jet);

% surf(x(1):precission:x(2),y(1):precission:y(2),a,'CDataMapping','direct')
% 
% hold on
% image(x,y,a,'CDataMapping','direct');

contour3(x(1):precission:x(2),y(1):precission:y(2),a,[0,0],'LineWidth',2,'Color','black')

% image(x,y,a,'CDataMapping','scaled','AlphaData',0.5);

plot3(train1(1,:),train1(2,:),log1,'b.','color','c');
plot3(train2(1,:),train2(2,:),log2,'b.','color','r');

colorbar();

xlim(x);
ylim(y);
zlim(z);

view(145,35);

% % % % % % % % % % % % errorcount = 0;
% % % % % % % % % % % % 
% % % % % % % % % % % % runns = 10000;
% % % % % % % % % % % % 
% % % % % % % % % % % % for i=1:runns/2
% % % % % % % % % % % % % testword = 1;
% % % % % % % % % % % % % hmm_example_generate_testseq
% % % % % % % % % % % % x2=simulate_test_sample(1);
% % % % % % % % % % % % [class,tmp]=bayes_classifier(param1,param2,[0.5,0.5],x2);
% % % % % % % % % % % % if(class == 2)
% % % % % % % % % % % % errorcount = errorcount +1;
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % % testword = 2;
% % % % % % % % % % % % % hmm_example_generate_testseq
% % % % % % % % % % % % x2=simulate_test_sample(2);
% % % % % % % % % % % % [class,tmp]=bayes_classifier(param1,param2,[0.5,0.5],x2);
% % % % % % % % % % % % if(class == 1)
% % % % % % % % % % % % errorcount = errorcount +1;
% % % % % % % % % % % % end
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % errorrate = errorcount/(runns)