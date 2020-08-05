function assert(expression, message)

% function assert(expression, message)
%
% Inputs:
%	expression = value of assertion (numeric)
%	message = an optional string
% Outputs:
% Description:
%	Check assertion value for debugging purposes.  If expression is false,
%	report error and die.

%-------------------------------------------------------------------------------
%
% File: assert.m
% Description:
% Author: Paul P. Audi, Engineer BBN
%
% $Date: 2002/06/25 16:46:23 $
% $Revision: 1.1 $
% $Id: assert.m,v 1.1 2002/06/25 16:46:23 class Exp $
%
% Revision History:
%
%-------------------------------------------------------------------------------

if (~expression),
  if (nargin == 1 | isempty(message)),
    error('no message');
  else,
    error(message);
  end;
end;
