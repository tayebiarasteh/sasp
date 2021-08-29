%Created on June 2016.

%https://github.com/tayebiarasteh/
%%
function [x,n] = impseq(n0,n1,n2)
% Impulse signal
% Generates x(n) = delta(n-n0); n1 <= n <= n2
% ----------------------------------------------
% [x,n] = impseq(n0,n1,n2)
%
n = n1:n2;
x = (n-n0) == 0;
