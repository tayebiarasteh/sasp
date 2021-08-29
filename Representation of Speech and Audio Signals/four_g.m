%Created on June 10, 2020.

%@author: Soroosh Tayebi Arasteh <soroosh.arasteh@fau.de>
%https://github.com/tayebiarasteh/

%% 2.4.g) Analysis in the Frequency Domain
clear
clc

v2 = randn(1,10000);
v4 = 0.5*v2;

V4 = (abs(fft(v4)).^2)/10000;

plot(V4)
