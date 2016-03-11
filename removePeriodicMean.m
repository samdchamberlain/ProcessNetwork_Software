function [anomalyMat]=removePeriodicMean(signalMat,period,sliderWidth,NoDataCode)

%Removes the mean at every timestep from a periodic signal.  For example,
%For several years of daily data, feed in a signal with daily resolution
%and give the "period" parameter as 365.  This will give you the daily
%anomaly.  The slider width is applied to limit the mean to a certain
%number of periods after the data point - making each value the anomaly
%from the mean of the next "sliderWidth" values at the same position in
%subsequent periods.

%signalMat is an n-column matrix where n variables are represented as a
%timeseries- columns are variables, rows are timeseries records

%period is the length of the repeating pattern in the data- if using
%FluxNet data with 30 min resolution, the period is 48 or one day

%sliderWidth is the number of periods to base the anomaly on - for a five
%day moving anomaly, set sliderWidth to 5.

%NoDataCode is a data value that will be ignored. 

signalMat(signalMat == NoDataCode) = NaN;

% Initialize
sz = size(signalMat);
aveMat = NaN(sz);

% Compute moving-periodic average
for i = 1:sz(1)

    % Get averaging indices
    ai = i-period*floor(sliderWidth/2):period:i+period*floor(sliderWidth/2);
    
    % Handle edges
    if min(ai) < 1
        ai = ai+floor(1-min(ai/period))*period;
    elseif max(ai) > sz(1)
        ai = ai-(ceil((max(ai)-sz(1))/period)*period);
    end
    
    % Average
    aveMat(i,:) = nanmean(signalMat(ai,:),1);
    
end

% Compute anomaly
anomalyMat = signalMat-aveMat;
anomalyMat(isnan(anomalyMat)) = NoDataCode;

