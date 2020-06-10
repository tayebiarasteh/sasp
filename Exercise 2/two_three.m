%Created on June 10, 2020.

%@author: Soroosh Tayebi Arasteh <soroosh.arasteh@fau.de>
%https://github.com/starasteh/

%% 2.2
clear
clc

[x1, fs1]=audioread('female_german.wav');
[x2, fs2]=audioread('male_german.wav');

% averaging over the channels, stereo to mono
x1 = mean(x1, 2);
x2 = mean(x2, 2);

% plotting based on seconds
figure
plot(0:1/fs1:(length(x1)-1)/fs1, x1)
title('female german')
xlabel('sec')

figure
plot(0:1/fs2:(length(x2)-1)/fs2, x2)
title('male german')
xlabel('sec')
grid

%% 2.3

% sending the sound to the speaker
sound(x1, fs1)
sound(x2, fs2)

% resampling
x1_mod = resample(x1, 8e3, fs1);
sound(x1_mod, 8e3)
sound(x1, fs1)
