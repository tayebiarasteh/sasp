function parm=ann_hmm(parm,var_fac,A_fac);
% function parm=ann_hmm(parm,var_fac,A_fac);

N=parm.N;
format compact

parm.Pi = parm.Pi + (A_fac-1) * 1/N;
parm.Pi = parm.Pi / sum(parm.Pi);
parm.A=parm.A + (A_fac-1)/N;
parm.A=parm.A ./ repmat(sum(parm.A,2),1,N);
for i=1:N,
      nmode = length(parm.pdf(i).modes);
      for j=1:nmode,
        parm.pdf(i).modes(j).cholesky_covar = ...
             parm.pdf(i).modes(j).cholesky_covar * var_fac;
        parm.pdf(i).modes(j).weight = parm.pdf(i).modes(j).weight + .001;
      end;
end;


