function decision_value = calculate_decision(limX,limY,resolution,gmm1,gmm2,priors)
%CALCULATE_DECISSION Summary of this function goes here
%   Detailed explanation goes here
x = limX(1):resolution:limX(2);
y = limY(1):resolution:limY(2);

coord(1,:) = reshape((ones(length(y),1)*x)',1,[]);
coord(2,:) = reshape((ones(length(x),1)*y),1,[]);

logPosterior1 = reshape(log10(pdf(gmm1,coord')),length(x),length(y))';
logPosterior2 = reshape(log10(pdf(gmm2,coord')),length(x),length(y))';

decision_value = logPosterior1-logPosterior2+log10(priors(1))-log10(priors(2));

end

