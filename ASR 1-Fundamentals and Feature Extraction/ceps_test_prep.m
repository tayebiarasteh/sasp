% ===============================================
%
% ILLUSTRATION OF CEPSTRAL ANALYSIS AND LIFTERING
%
% Herbert Buchner, Jan 2003
%
% ===============================================
echo off

% load speech signal
%[y,fs] = wavread('tapestry.wav'); % female

[y,fs] = wavread('male_german.wav'); % male
y=resample(y,16000,fs);
y = y(9500:end,1);
fs=16000;

% how many samples correspond to 20ms block at sampling rate fs ?
blocklen = fs * 20e-3;

fftlen = 4096; % long FFT for better visualization

for i = 1:floor(length(y)/blocklen)

  % === 1. ORIGINAL SIGNAL BLOCK ===
  block = y((i-1)*blocklen+1:i*blocklen);

  block_w = block.*hamming(blocklen);

  subplot(411)
  plot(block);
  title('ORIGINAL SIGNAL BLOCK')

  % === 2. LOG POWER SPECTRAL DENSITY ===
  subplot(412)
  block_psd=20*log10(abs(fft(block_w,fftlen)));
  plot(block_psd(1:fftlen/2))
  title('LOG POWER SPECTRAL DENSITY')

  % === 3. REAL CEPSTRUM ===
  subplot(413)
%  block_ceps= ...  ; % == FILL IN THE FORMULA. Use fft of length fftlen.
%  plot(block_ceps(1:blocklen))
  title('REAL CEPSTRUM')

  % === 4. LIFTERED LOG POWER SPECTRAL DENSITY USING COMPLEX CEPSTRUM ===
  subplot(414)
  % complex cepstrum
%  block_cceps= ...  ; % == FILL IN THE FORMULA. Use fft of length fftlen.
%  liftered_block_cceps= ...  ; % == FILL IN THE FORMULA.
%  liftered_block_psd=20*log10( ... ); % == FILL IN THE FORMULA. Use fft of length fftlen.
%  plot(liftered_block_psd(1:fftlen/2))
  title('LIFTERED LOG POWER SPECTRAL DENSITY USING COMPLEX CEPSTRUM')

  pause

end

