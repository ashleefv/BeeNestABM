function params = assignCohortParameters(COHORTparams,params,tags)
    params(tags == 0) = COHORTparams(2); % untreated
    params(1) = COHORTparams(1); % queen
    params(tags == 1) = COHORTparams(3); % sucrose control
    params(tags == 2) = COHORTparams(4); % low dose
    params(tags == 3) = COHORTparams(5); % high dose
end
