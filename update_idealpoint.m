function idealpoint = update_idealpoint(idealpoint,population)
fitnessvec = [population(:).Fit]';
idealpoint = min([idealpoint;fitnessvec],[],1);

end