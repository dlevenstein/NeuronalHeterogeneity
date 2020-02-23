function plogt = LogGamma(log10rate,CV,weight,tau)
    %Parameters: 1) log10 rate  2) log10CV . 3) weights
    %Tau must be in first dimension. other parameters can be as many
    %dimensions as you want as long as first is singleton?
    k = 1./(CV); %alpha
    %lambda = exp(lambkweit(1:end/3));  %beta
    lambda = (10.^log10rate).*k;
    
    if length(k)==1 && length(lambda)==1 && length(weight)>1
        %If multiple weights, make a big matrix... (shared mode)
        k = k.*ones(size(weight));
        lambda = lambda.*ones(size(weight));
    end
    
    %If multiple taus and parameters, make a big matrix...
    k = repmat(k,length(tau),1);
    lambda = repmat(lambda,length(tau),1);
    weight = repmat(weight,length(tau),1);

    plogt = weight.*(lambda.^(k).*exp(k.*tau)) ./ (gamma(k).*exp(lambda.*exp(tau)));
    %plogt(isnan(plogt)) = 0;
end

