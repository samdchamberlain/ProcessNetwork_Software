function [Surrogates] = createSurrogates(opts,Data,nsur)
% Create the surrogate data for statistical testing
% 
% ----------- Inputs -----------
% opts...
%
% ---------- Outputs -----------
% Surrogates = the surrogate data
%
% ------------------------------

Surrogates = NaN(size(Data,1),size(Data,2),nsur);

if opts.SurrogateMethod == 2

    % Randomly shuffle data (leaving NaNs where they are)
    for i = 1:size(Data,2)
        ni = ~isnan(Data(:,i));
        for ti = 1:nsur
            Surrogates(ni,i,ti) = randsample(Data(ni,i),sum(ni));
        end
    end
elseif opts.SurrogateMethod == 3
    
    % Create Iterated Amplitude Adjusted Fourier Transform surrogates
    if isempty(Data(Data == opts.NoDataCode))
        for i = 1:size(Data,2)
            Surrogates(:,i,:) = IAAFTsur(Data(:,i),nsur);
        end
    else
        % IAAFT surrogates requires gap-free data
        logwrite('Warning: Surrogates set to NaN. Use of IAAFT method requires gap-free data.',1);
    end
elseif opts.SurrogateMethod == 4
    
    % Create random walks
    if isempty(Data(Data == opts.NoDataCode))
        for ti = 1:nsur
            Surrogates(:,:,ti) = cumsum(randn(size(Data,1), size(Data,2)));
        end    
    else
        % Random walks used for wavelets, which requires gap-free data
        logwrite('Warning: Surrogates set to NaN. Use gap-free data for random walks with wavelets.',1);
    end
   
end
