% Exercise 6: Human Machine Interfaces
% SLOC, BF and AEC
% Lutz Marquardt & Edwin Mabande 
% LMS Audio, 01.10

clear

% Initialisation
d = 0.04;                                   % inter microphone distance
M = 5;                                      % number of microphones
distance_outer_mics = d*(M-1);              % array length
c = 342;                                    % speed of sound

% loading 5 channel wav files
[mic_D, fs] = audioread('Desired_speech.wav'); % microphone signals with only the desired speaker
mic_I = audioread('Interferer_speech.wav');    % microphone signals with only the interfering speaker

sig_len = length(mic_D);                    % (30 seconds)
tmax = ceil(d*(M-1)/c*fs);                  % Maximum searched TDOA (in samples)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        SLOC                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Estimate TDOA using CC function')

% Compute the crosscorrelation between microphone signals (mic 0 and mic4),
% for desired and interfering speaker (use maximum lag: tmax). Determine the
% respective delays in samples. 

% Crosscorrelations:
D_rxjxi = ...
I_rxjxi = ...

% TDOAs: 
[D_val, D_tij_estim] = ...         % D_tij_estim should be between -tdoa_max and +tdoa_max
[I_val, I_tij_estim] = ...         % I_tij_estim should be between -tdoa_max and +tdoa_max

% Plot 
D_tij_true = 8;
figure
subplot(1,2,1)
plot([-tmax:tmax],D_rxjxi);
hold on
plot(D_tij_estim,D_val,'or','Markersize',20,'Linewidth',4)
title({['True TDOA: ' num2str(D_tij_true) 'samples'];['Estimated TDOA: ' num2str(D_tij_estim) 'samples']});
xlabel('samples')
ylabel('CC function')
grid on

I_tij_true = -8;
subplot(1,2,2)
plot([-tmax:tmax],I_rxjxi);
hold on
plot(I_tij_estim,I_val,'or','Markersize',20,'Linewidth',4)
title({['True TDOA: ' num2str(I_tij_true) 'samples'];['Estimated TDOA: ' num2str(I_tij_estim) 'samples']});
xlabel('samples')
ylabel('CC function')
grid on


% Compute the TDOAs (in seconds) and the corresponding DOAs (in degrees)

TDOA = ...
theta_D = ...

TDOA = ...
theta_I = ...

% Compute the delay in samples between adjacent microphones 

dly_adj_mikes = ...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             BEAMFORMER                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Delay-and-Sum Beamformer');

% creating composite microphone signals: contain desired and interferer speech
mic_C = mic_D + mic_I;

a_m=1/M; % uniform weighting factor

% Implement a DSB steering into the direction of the desired speaker, using
% the DOA computed by the source localization

% DSB output for desired signal, interferer and composite signal:
bf_out_desired = ...
bf_out_interferer = ...
bf_out_composite = ...

% % Signal Playback
% % Desired
% sound(mic_D(:,1),fs);                         % reference microphone
% sound(bf_out_desired,fs);                     % beamformer output
% % Interferer
% sound(mic_I(:,1),fs);                         % reference microphone
% sound(bf_out_interferer,fs);                  % beamformer output
% % Combined
% sound(mic_C(:,1),fs);                         % reference microphone
% sound(bf_out_composite,fs);                     % beamformer output


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   AEC                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Acoustic Echo Cancellation')

[y,fs_down] = audioread('bf_aec.wav');
speaker = audioread('Speaker.wav');

% algorithmic parameters
N_hat = 2000;   % filter length 
mu = 0.5;       % stepsize
delta = 0.001;  % regularization parameter

% initialization
L=length(y);
e=zeros(L,1);               % error signal
xi=zeros(L,1);              % xi=0 if double talk is present, xi=1 if not
xi=((1:L)<double_talk_start)+((1:L)>double_talk_end); % perfect doubletalk detector

h_hat=zeros(N_hat,1);       % AEC-filter
xvec=zeros(N_hat,1);        % Reference signal input to AEC-filter
x=resample(speaker,1,4);    % 12 kHz speaker signal

% adaptation
for i=1:L
    xvec = ...
    y_hat = ...
    e(i) = ...              % error signal computation
    h_hat = ...             % update equation
    % disp(i/L)             % display computing progress
end


figure
% ...                       % plot error signal e
grid on

% % Signal Playback
% sound(e,fs_down);


