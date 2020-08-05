% function i = select_discrete_rv(Pi)
%
% Inputs:
%	Pi = vector of probabilities
% Outputs:
%	i = discrete random variable
% Description:
%	This function selects a discrete random variable in 1:length(Pi)
%	based on the probabilities in Pi.
%	The probabilities in Pi should be non-negative and sum to 1.0

%-------------------------------------------------------------------------------
%
% File: select_discrete_rv.m
% Description:
% Author: Dr. Paul M. Baggenstoss, NUWC, Newport, RI, Code 21
%
% $Date: 2002/06/25 16:45:36 $
% $Revision: 1.1 $
% $Id: select_discrete_rv.m,v 1.1 2002/06/25 16:45:36 class Exp $
%
% Revision History:
%
%-------------------------------------------------------------------------------

function i = select_discrete_rv(Pi)

assert(all(Pi >= 0.0));

n_i = length(Pi);
val = rand;
sum = 0;
for i=1:n_i,
    sum = sum + Pi(i);
    if(val < sum), return; end;
end;
return

