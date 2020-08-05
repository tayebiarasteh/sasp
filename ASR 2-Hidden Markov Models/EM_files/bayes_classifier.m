function [class,logLikeliehood,varargout] = bayes_classifier(gmm1,gmm2,prior,data)
%BAYES_CLASSIFIER Summary of this function goes here
%   Detailed explanation goes here
    logPosterior1 = log10(pdf(gmm1,data'));
    logPosterior2 = log10(pdf(gmm2,data'));
    
    logLikeliehood1 = ones(1,length(logPosterior1))*logPosterior1+log10(prior(1));
    logLikeliehood2 = ones(1,length(logPosterior2))*logPosterior2+log10(prior(2));
    
    logLikeliehood = [logLikeliehood1,logLikeliehood2];
    
    [tmp,class]=max(logLikeliehood);
    
    %for classifier plot
    
    if(nargout == 3)
    
    logPosterior1 = log10(pdf(gmm1,data'));
    logPosterior2 = log10(pdf(gmm2,data'));

    varargout{1} = logPosterior1-logPosterior2+log10(prior(1))-log10(prior(2));
    end

end

