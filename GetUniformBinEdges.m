function [binEdges,minEdge,maxEdge]=GetUniformBinEdges(sampleMat,nBinMat,pctlRange,NoDataCode)

[~,nSignals]=size(sampleMat);
binMax=max(nBinMat);
binEdges=NaN(nSignals,binMax);
minEdge=NaN(nSignals,1);
maxEdge=NaN(nSignals,1);

%make nodata entries NaN, because min and max function ignores them...
sampleMat(sampleMat == NoDataCode) = NaN;

%compute the bin edges using fractions of the min and max
for s=1:nSignals
        % Use uniform bin width
        minEdge(s)=prctile(sampleMat(:,s),pctlRange(1));
        maxEdge(s)=prctile(sampleMat(:,s),pctlRange(end));
        E = linspace(minEdge(s),maxEdge(s),nBinMat(s)+1); % all edges, including start
        binEdges(s,1:nBinMat(s))= E(2:end);
end