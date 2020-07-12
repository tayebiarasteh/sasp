% 3.3.1. Masking in the Frequency Domain
clear
clc

%% a
% fs = 8000;
% dt = 1/fs;
% duration = 1;
% t = 0:dt:duration;
% F = 1000; % fundamental frequency
% masker = sin(2*pi*F*t);
% 
% stem(t,masker)
% sound(masker)
%% b

t = 0:1/1e3:2;
y = chirp(t,0,1,250);

plot(y)