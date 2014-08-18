function threshold = thresholdDetermination(sampSizes)
%thresholdDetermination calculates threshold for numeric vector
% 
% 
    
% get 90 percent quantile
    upperQuantile = quantile(sampSizes, 0.9);
    
    % scale quantile down to get threshold
    threshold = 0.9 * upperQuantile;
    
end