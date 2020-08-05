% =========================================
%
% MFCC feature vectors for different vocals
%
% Herbert Buchner, Jan 2003
%
% =========================================


load vocals_herbert
%load vocals_robert
figure
grid on
hold on

[ceps_a,freqresp,fb,fbrecon,freqrecon] = mfcc(ya,16000,100);
plot(ceps_a(2,:),ceps_a(5,:),'b.')
pause

[ceps_i,freqresp,fb,fbrecon,freqrecon] = mfcc(yi,16000,100);
plot(ceps_i(2,:),ceps_i(5,:),'g.')
pause

[ceps_u,freqresp,fb,fbrecon,freqrecon] = mfcc(yu,16000,100);
plot(ceps_u(2,:),ceps_u(5,:),'r.')
pause

[ceps_e,freqresp,fb,fbrecon,freqrecon] = mfcc(ye,16000,100);
plot(ceps_e(2,:),ceps_e(5,:),'c.')
pause

[ceps_o,freqresp,fb,fbrecon,freqrecon] = mfcc(yo,16000,100);
plot(ceps_o(2,:),ceps_o(5,:),'m.')
pause

hold off

