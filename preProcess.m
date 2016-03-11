function [Data] = preProcess(opts,Data)
% Run the preprocessing options, currently includes data trimming and
% transformation
% 
% ----------- Inputs -----------
% opts...
%
% ---------- Outputs -----------
% rawData = the preprocessed data
%
% ------------------------------

% Turn NoDataCode into NaN
Data(Data == opts.NoDataCode) = NaN;

% Trim rows where there is missing data
if opts.trimTheData
    Data(sum(isnan(Data),2) > 0,:) = NaN;
end

% Apply transformation
if opts.transformation == 1
    
    % Apply anomaly filter 
    Data = removePeriodicMean(Data,opts.anomalyPeriodInData,opts.anomalyMovingAveragePeriodNumber,NaN);

elseif opts.transformation == 2
    
    % Apply wavelet transform
    if sum(isnan(Data(:))) == 0
        Data = waveletTransform(Data,opts.waveN,opts.waveName,opts.waveDorS,opts.parallelWorkers);
    else
       logwrite('Warning: Wavelet transform not applied. Gap-free data required.',1);
    end

end
