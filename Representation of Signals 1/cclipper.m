%Created on June 6, 2020.

%@author: Soroosh Tayebi Arasteh <soroosh.arasteh@fau.de>
%https://github.com/starasteh/
%%
function y = cclipper(x,eta)

for k = 1:length(x)
	if x(k) > eta
		y(k) = x(k);
	elseif x(k) < -eta
		y(k) = x(k);
	else
		y(k) = 0;
	end
end