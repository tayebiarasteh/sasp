%Created on June 10, 2020.

%@author: Soroosh Tayebi Arasteh <soroosh.arasteh@fau.de>
%https://github.com/tayebiarasteh/

%% 2.4.a) Long-Term and Short-Time Analysis
clear
clc

[x1, fs1]=audioread('female_german.wav');
[x2, fs2]=audioread('male_german.wav');
% averaging over the channels
x1 = mean(x1, 2);
x2 = mean(x2, 2);

% long-term pdf
subplot(2,1,1)
hist(x2, 100)
title('entire male signal')

% short-term pdf
subplot(2,1,2)
hist(x2(1000:1000+0.025*fs2), 100)
title('25ms male signal')

%% b)

[x3, fs3]=audioread('castanets.wav');
[x4, fs4]=audioread('orchestra.wav');
% averaging over the channels, stereo to mono
x3 = mean(x3, 2);
x4 = mean(x4, 2);

subplot(2,1,1)
hist(x3, 100)
title('castanets')

subplot(2,1,2)
hist(x4, 100)
title('orchestra')

%% c) Analysis in the Time Domain

v1 = randn(1,1000);
v2 = randn(1,10000);

% ACF
Rv1v1_biased = xcorr(v1,'biased');
Rv1v1_unbiased = xcorr(v1,'unbiased');
Rv2v2_biased = xcorr(v2,'biased');
Rv2v2_unbiased = xcorr(v2,'unbiased');

subplot(2,2,1)
plot(Rv1v1_biased)
title('Rv1v1 biased')
subplot(2,2,2)
plot(Rv1v1_unbiased)
title('Rv1v1 unbiased')
subplot(2,2,3)
plot(Rv2v2_biased)
title('Rv2v2 biased')
subplot(2,2,4)
plot(Rv2v2_unbiased)
title('Rv2v2 unbiased')

%% e) adding a mean of one
v3 = v2 + 1;
Rv3v3_biased = xcorr(v3,'biased');
Rv3v3_unbiased = xcorr(v3,'unbiased');

subplot(2,1,1)
plot(Rv3v3_biased)
title('Rv3v3 biased')
subplot(2,1,2)
plot(Rv3v3_unbiased)
title('Rv3v3 unbiased')

%% f)

figure
Rx1x1_biased = xcorr(x1(1:fs1*6),'biased');
subplot(2,2,1)
plot(Rx1x1_biased)
title('Rx1x1 biased')
Rx2x2_biased = xcorr(x2(1:fs2*6),'biased');
subplot(2,2,2)
plot(Rx2x2_biased)
title('Rx2x2 biased')
Rx3x3_biased = xcorr(x3(1:fs3*6),'biased');
subplot(2,2,3)
plot(Rx3x3_biased)
title('Rx3x3 biased')
Rx4x4_biased = xcorr(x4(1:fs4*6),'biased');
subplot(2,2,4)
plot(Rx4x4_biased)
title('Rx4x4 biased')

figure
Rx1x1_biased = xcorr(x1(100000:100000+fs1*0.02),'unbiased');
subplot(2,2,1)
plot(Rx1x1_biased)
title('Rx1x1 unbiased')
Rx2x2_biased = xcorr(x2(100000:100000+fs2*0.02),'unbiased');
subplot(2,2,2)
plot(Rx2x2_biased)
title('Rx2x2 unbiased')
Rx3x3_biased = xcorr(x3(100000:100000+fs3*0.02),'unbiased');
subplot(2,2,3)
plot(Rx3x3_biased)
title('Rx3x3 unbiased')
Rx4x4_biased = xcorr(x4(100000:100000+fs4*0.02),'unbiased');
subplot(2,2,4)
plot(Rx4x4_biased)
title('Rx4x4 unbiased')
