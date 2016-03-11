function [binEdgesGlobal,minEdgeGlobal,maxEdgeGlobal]=GetEvenBinEdgesGlobal(nBinVect,minEdgeLocal,maxEdgeLocal)
nVars = size(maxEdgeLocal,1);

minEdgeGlobal = nanmin(minEdgeLocal,[],2);
maxEdgeGlobal = nanmax(maxEdgeLocal,[],2);
binEdgesGlobal = NaN(nVars,max(nBinVect));
for i = 1:nVars
    E = linspace(minEdgeGlobal(i),maxEdgeGlobal(i),nBinVect(i)+1); % all edges, including start
    binEdgesGlobal(i,1:nBinVect(i))= E(2:end);
end
