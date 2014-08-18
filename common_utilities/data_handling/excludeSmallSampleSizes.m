function [remainingData, excluded] = excludeSmallSampleSizes(dataStruct)
% EXCLUDESMALLSAMPLESIZES - exclude assets with too less observations
%   
% Args:
%  dataStruct               1 x nAss structure array as returned
%                           by hist_stock_data
% Output:
%  remainingData            1 x someAss structure array as would
%                           be returned by hist_stock_data, if
%                           assets with too less observations
%                           had not been included into the
%                           download
%
%  excluded                 structure with information about excluded
%                           assets, with sample sizes and relative sample
%                           sizes

    
% get number of asset
    nAss = length(dataStruct);
    
    % get sample sizes
    sampSizes = ones(1,nAss);
    for ii=1:nAss
        sampSizes(ii) = length(dataStruct(ii).Date);
    end
    
    % get threshold in relative terms to upper quantile 
    threshold = thresholdDetermination(sampSizes);
    upperQuantile = threshold / 0.9;
    
    % find assets with too few observations
    shouldAssetBeExcluded = sampSizes < threshold;
    assetsToExclude = find(shouldAssetBeExcluded);
    nExcluded = sum(shouldAssetBeExcluded);
    
    % we need for loop to access same field for multiple components
    % in a structure array
    excludedAssetNames = cell(1, nExcluded);
    excludedSampleSizes = zeros(1, nExcluded);
    excludedSampleSizesRelative = zeros(1, nExcluded);
    
    % get information on excluded variables
    for ii=1:nExcluded
        % index of current asset that will be excluded
        ind = assetsToExclude(ii);
        
        % get information about asset
        excludedSampleSizes(ii) = sampSizes(ind);
        excludedSampleSizesRelative(ii) = ...
            sampSizes(ind) / upperQuantile;
        excludedAssetNames{ii} = dataStruct(ind).Ticker;
    end
    
    % summarize information in structure 
    excluded = struct('Assets',  {excludedAssetNames}, ...
                      'Sample_sizes', excludedSampleSizes, ...
                      'Relative_sample_sizes', ...
                      excludedSampleSizesRelative);
    
    % get remaining assets with sufficient data
    dataStruct(assetsToExclude) = [];
    
    remainingData = dataStruct;
    
end

