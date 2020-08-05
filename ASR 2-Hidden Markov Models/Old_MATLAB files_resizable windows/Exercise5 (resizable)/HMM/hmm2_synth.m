function [feats,states] = hmm2_synth(hmm_parm, model_param, nsamp);
% function [feats,states] = hmm2_synth(hmm_parm, model_param, nsamp);

if( nargin < 3 ),
  error('Usage: [feats,states] = hmm2_synth(hmm_parm, model_param, nsamp);');
end;

[feats,states]=hmm2_synth_mex(hmm_parm,nsamp);
feats=feats';

return

