echo on
% ===========================================================
%
% ILLUSTRATION OF MEL-FREQUENCY CEPSTRUM COEFFICIENTS FOR ASR
%
% Herbert Buchner, Jan 2003
%
% ===========================================================
%
figure

% Load utterance from TIMIT database
% 'A huge tapestry hung in her hallway'
% (TRAIN/DR5/FCDR1/SX106).
[tap,fs]=wavread('tapestry.wav');
fs
length(tap)

subplot(211)
plot(tap)
axis([1 length(tap) min(tap) max(tap)]);

% Calculate Mel Frequency Cepstrum Coefficients (MFCC) and some
% intermediate results that can be used for analyzing the procedure.
% Sampling rate 16kHz, pictures are sampled at 100Hz. There are 312 frames.
[ceps,freqresp,fb,fbrecon,freqrecon] = mfcc(tap,16000,100);

pause


% As an intermediate result, the (uncompressed) DFT spectrogram (freqresp) is
% first calculated. (It has been flipped so that high frequencies are at the top.)
subplot(212)
imagesc(flipud(freqresp));
colormap(1-gray);
axis([1 size(freqresp,2) 1 size(freqresp,1)]);

pause


% After combining several DFT channels into a single mel-scale channel, the
% result is the filter bank output. (The fb output of the mfcc command includes
% the log10 calculation.)
imagesc(flipud(fb));
axis([1 size(fb,2) 1 size(fb,1)]);

pause


% Finally, the conversion into the cepstral domain uses the discrete cosine
% transform to reduce the dimensionality of the output.
% This leads to the 12 MFCC coefficients C1 - C12.
imagesc(ceps(2:end,:))
axis([1 size(ceps,2) 1 12]);

pause

% Finally, this set of MFCC coefficients is then extended by the
% coefficients C0 which are a functin of the input signal power.
% Since the wave form in our work is normalized to be between -1 and +1, the
% C0 coefficients are all negative. The other coefficients, C1 - C12, are
% generally zero mean.
plot(ceps(1,:))
pause

imagesc(ceps);

pause


% ===========================================================================


% We can invert the cosine transform to get back
% into the filter bank domain and see
% how much information we have lost. The recon
% output is again flipped so that the high frequency channel is at the top.
% Note the result is much smoother than the original filter bank output.
subplot(211)
imagesc(flipud(fb));
subplot(212)
imagesc(flipud(fbrecon));

pause


% The original spectrogram is reproduced by resampling the result above.
% Note the fine lines indicating the pitch have been lost by the MFCC
% process. This is a useful feature for speech recognition
% systems since they want a representation that is blind to pitch changes.
subplot(211)
imagesc(flipud(freqresp));
subplot(212)
imagesc(flipud(freqrecon));

%pause

%figure
%plot(ceps(2,:),ceps(5,:))
%%plot3(ceps(2,:),ceps(5,:),ceps(7,:))




