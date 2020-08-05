% intermediate functions

%% 1.3. center clipping
clear
clc

k = -32:0.2:32;
v12 = sin(2*pi*k/32);
y12 = cclipper(v12,0.1);
subplot(2,1,1)
plot(k, v12)
ylabel('v12[k]')
grid on
subplot(2,1,2)
plot(k, y12)
ylabel('y12[k]')
grid on
%% 1.4.1.b
clear
clc

x1 = rand(1,1000);
x2 = rand(1, 10000);
subplot(2,2,1)
plot(x1);
subplot(2,2,2)
plot(x2);
mx1 = mean(x1);
mx2 = mean(x2);
varx1 = var(x1);
varx2 = var(x2);
subplot(2,2,3)
hist(x1, 100);
subplot(2,2,4)
hist(x2, 100);
%% 1.4.1.c
clear
clc

x3 = randn(1,1000);
x4 = randn(1, 10000);
subplot(2,2,1)
plot(x3);
subplot(2,2,2)
plot(x4);
mx1 = mean(x1);
mx2 = mean(x2);
varx1 = var(x3);
varx2 = var(x4);
subplot(2,2,3)
hist(x3, 100);
subplot(2,2,4)
hist(x4, 100);
%% 1.4.3.b
clear
clc

b1 = [ 0.5 0.5 ];
c1 = [ 1 0 ];
b2 = [ 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 ];
c2 = [ 1 0 0 0 0 0 0 0 0 0 ];
x2 = rand(1, 10000);
x4 = randn(1, 10000);

% filtering
y21 = filter(b1, c1, x2);
y22 = filter(b2, c2, x2);
y41 = filter(b1, c1, x4);
y42 = filter(b2, c2, x4);
subplot(2,2,1)
hist(y21, 100);
title('y21, uniform')
subplot(2,2,2)
hist(y22, 100);
title('y22, uniform')
subplot(2,2,3)
hist(y41, 100);
title('y41, Gaussian')
subplot(2,2,4)
hist(y42, 100);
title('y42, Gaussian')

% ACF
subplot(2,3,1)
x2_ACF = xcorr(x2);
plot(x2_ACF)
title('x2 ACF, uniform')
subplot(2,3,2)
y21_ACF = xcorr(y21);
plot(y21_ACF)
title('y21 ACF, uniform')
subplot(2,3,3)
y22_ACF = xcorr(y22);
plot(y22_ACF)
title('y22 ACF, uniform')
subplot(2,3,4)
x4_ACF = xcorr(x4);
plot(x4_ACF)
title('x4 ACF, Gaussian')
subplot(2,3,5)
y41_ACF = xcorr(y41);
plot(y41_ACF)
title('y41 ACF, Gaussian')
subplot(2,3,6)
y42_ACF = xcorr(y42);
plot(y42_ACF)
title('y42 ACF, Gaussian')

% PSD
subplot(2,3,1)
x2_PSD = pwelch(x2);
plot(x2_PSD)
title('x2 PSD, uniform')
subplot(2,3,2)
y21_PSD = pwelch(y21);
plot(y21_PSD)
title('y21 PSD, uniform')
subplot(2,3,3)
y22_PSD = pwelch(y22);
plot(y22_PSD)
title('y22 PSD, uniform')
subplot(2,3,4)
x4_PSD = pwelch(x4);
plot(x4_PSD)
title('x4 PSD, Gaussian')
subplot(2,3,5)
y41_PSD = pwelch(y41);
plot(y41_PSD)
title('y41 PSD, Gaussian')
subplot(2,3,6)
y42_PSD = pwelch(y42);
plot(y42_PSD)
title('y42 PSD, Gaussian')
