% Basic Linear Predictive Coding (forward),
% Herbert Buchner, Jan. 2002
% ==========================================================

% ------------
% input signal
% ------------

s=audioread('male.wav');  % sampling rate 8kHz



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

for m=2:no_iterations

  % windowed frame (including samples from past frame)
  s_p=hamming(L1+L2).*s( L1.*(m-1)-L2 : L1.*m-1 );

  a=real(lpc(s_p,N));    % N+1 by 1 coefficient vector
  A=[A; a];

  % convolution
  for k=0:L1-1
    err(L1.*(m-1)+k)=a*s(L1.*(m-1)+k:-1:L1.*(m-1)-N+k);
  end

end


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









