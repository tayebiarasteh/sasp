% Basic Linear Predictive Coding (forward),
% Christian Huemmer, 2015
% ==========================================================

% ------------
% input signal
% ------------

s=audioread('male.wav'); s=s(1:32000); %four seconds of speech 
fs=8000; % sampling rate 8kHz


% ----------
% LP encoder
% ----------

L1=160;  % block length
L2=80;   % overlap (previous block)
N=10;    % order of predictor

A=[];    % matrix for predictor coefficients
err=zeros(size(s,1),1);
y=zeros(size(s,1),1);
no_iterations =floor(size(s,1)/L1)-1;     % divided by L1 due to block processing
% Optionally: Visualize time-frequency dependent short-term PSD estimate
f=linspace(0,fs/2,L1+L2); % Vector with frequency bins
t=linspace(0,4,no_iterations); % Vector for the time index
PSD_s=zeros(no_iterations,length(f)); %PSD of input signal
PSD_err=zeros(no_iterations,length(f)); %PSD of output signal


for m=2:no_iterations

  % windowed frame (including samples from past frame)
  s_p=hamming(L1+L2).* s( L1.*(m-1)-L2 : L1.*m-1 );
  % time-dependent short-time PSD estimate
  PSD_s(m-1,:) = abs(fft(s_p)); 
  a=real(lpc(s_p,N));    % N+1 by 1 coefficient vector
  A=[A; a];

  % convolution
  for k=0:L1-1
    err(L1.*(m-1)+k)=a*s(L1.*(m-1)+k:-1:L1.*(m-1)-N+k);
  end
  PSD_err(m-1,:) = abs(fft( err(L1.*(m-1):L1.*(m-1)+L1-1),L1+L2)); 
end

%Plot time-frequency dependent short-time PSD
surf(t,f,PSD_s')
title('PSD of LPC input signal');
xlabel('time [sec]'); ylabel('frequency [Hz]'); zlabel('Ampl.');
figure;
surf(t,f,PSD_err')
title('PSD of LPC output signal');
xlabel('time [sec]'); ylabel('frequency [Hz]'); zlabel('Ampl.');
% -------------
% Quantization
% -------------

err2=quant(err,0.001);



% ----------
% LP decoder
% ----------
% Equation y(k)=err2(k)+a*y(k)

for m=2:no_iterations

  % convolution, synthesis filter
  for k=0:L1-1
    y(L1.*(m-1)+k)=err2(L1.*(m-1)+k)-A(m-1,2:N+1)*y(L1.*(m-1)+k-1:-1:L1.*(m-1)+k-N);
  end

end









