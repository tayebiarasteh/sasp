% function parm = init_gmix(data,n_modes,names,min_std,[rand_init]);
% Stand-alone version of gmix_init (does not require toolkit)
% data is dimensioned n_samp-by-n_feat
% names is a cellarray of strings

%  Date              Reason
%  ------------      --------------
%  Mar 28, 1997      Initial release
%  May 21, 1997      Handles non-integer nmode by rounding
%  Dec 10, 1999      Upgraded to MATLAB 5
%  Dec 20, 1999      took out normalization of min_std
%  Feb 27, 2001      Changed input arguments to be compatible with toolkit.
%  Apr 3, 2001       Disabled data_means, data_std and added check for ill-conditioning.

function parm = init_gmix(data,n_modes,names,min_std,rand_init);

if(nargin < 5), rand_init=1; end;

[n_feats,n_samples] = size(data);

parm = [];
for (i_feat = 1:n_feats),
  parm.features(i_feat).name = names{i_feat};
  parm.features(i_feat).min_std = min_std(i_feat);
end;

% Determine mean and std of input data
if( n_samples > 1 ),
  data_means = mean(data,2);
  data_std = std(data,0,2);
else,
  data_means = data;
  data_std = zeros(n_feats,1);
end;

% test for ill-conditioning
if(any( abs(data_means) > 1000*data_std) | ...
            max(data_std) > 100*min(data_std)),
       fprintf('***** GSUM_INIT: Red Alert!! Warning: Possible ill-conditioning: ');
       fprintf(' [data_means(:) data_std(:)]=\n');
       [data_means(:) data_std(:)]
       fprintf('Your features are shit. You should scale them or remove means *****\n');
       pause(1.0)
end;

if( n_samples > 1 ),
  xmin = min(data,[],2);
  xmax = max(data,[],2);
else,
  xmin = data;
  xmax = data;
end;

% The initial value of covariances
starting_std = max((xmax-xmin)/4,min_std(:));

if( rand_init),
  %--- Select starting means from a random set of data
  %--- Make sure there are none the same
  idx_has_duplicates = 1;
  n_modes_min = 1;
  while(idx_has_duplicates & n_modes >= n_modes_min),
    i_tries = 0;
    % try 1000 times
    while(idx_has_duplicates & i_tries < 1000),
      % generate n_modes random numbers in [1,n_samples]
      idx = 1 + floor(rand(1,n_modes) * n_samples);
      idx = sort(idx);
      % see if they are each unique
      if(length(idx) > 1),
	if( min( idx(2:n_modes) - idx(1:n_modes-1) ) > 0 ),
	  idx_has_duplicates = 0;
	end;
      else,
	idx_has_duplicates = 0;
      end;
      i_tries = i_tries + 1;
    end;
    % if this still fails after 1000 tries, reduce number of modes and keep trying
    if(idx_has_duplicates),
      n_modes = n_modes - 1;
      fprintf('gmix_init: reducing n_modes to %d,\n', n_modes);
      fprintf(' no. of samples (%d) too small compared to n_modes\n', n_samples);
    end;
  end;
  % number of modes has fallen below minimum, fail
  if(idx_has_duplicates),
    error('gmix_init: n_modes has been reduced to less than allowed, you need more data');
  end;
else,	% non-random mode mean initialization
  n_modes = min(n_modes,n_samples);
  idx = 1:n_modes;
end;

means = data(:,idx);

%  Create the (Cholesky of) covariances and store the means
wts = ones(1, n_modes) / n_modes;
for (i_mode = 1:n_modes),
  parm.modes(i_mode).cholesky_covar = diag( starting_std);
  %parm.modes(i_mode).cholesky_covar = eye(n_feats);
  parm.modes(i_mode).mean = means(:,i_mode);
  parm.modes(i_mode).weight = wts(i_mode);
end;

fprintf('init_gmix successful. Features:');
for i_feat = 1:n_feats,
   fprintf(' %s', parm.features(i_feat).name);
end;
fprintf('\n');

