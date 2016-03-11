
function [InormByDist,TnormByDist,SigThreshInormByDist,SigThreshTnormByDist,Ic,Tc,TvsIzero,SigThreshTvsIzero,RelI,RelT,HXtNormByDist,IvsIzero,SigThreshIvsIzero]=NormTheStats(nBinVect,I,T,SigThreshI,SigThreshT,meanShuffI,sigmaShuffI,meanShuffT,sigmaShuffT,HXt,HYf,lagVect)

% Get size
[nSignals,~,nLags,nFiles]=size(I);
nSLags = size(SigThreshI,3);

% Initialize
InormByDist=NaN(nSignals,nSignals,nLags,nFiles);
TnormByDist=NaN(nSignals,nSignals,nLags,nFiles);
SigThreshInormByDist=NaN(nSignals,nSignals,nSLags,nFiles);
SigThreshTnormByDist=NaN(nSignals,nSignals,nSLags,nFiles);
Ic=NaN(nSignals,nSignals,nLags,nFiles);
Tc=NaN(nSignals,nSignals,nLags,nFiles);
TvsIzero=NaN(nSignals,nSignals,nLags,nFiles);
SigThreshTvsIzero=NaN(nSignals,nSignals,nSLags,nFiles);
RelI=NaN(nSignals,nSignals,nLags,nFiles);
RelT=NaN(nSignals,nSignals,nLags,nFiles);
HXtNormByDist=NaN(nSignals,nSignals,nLags,nFiles);
IvsIzero=NaN(nSignals,nSignals,nLags,nFiles);
SigThreshIvsIzero=NaN(nSignals,nSignals,nSLags,nFiles);

l0i = find(lagVect == 0,1); % index of lag = 0

for f = 1:nFiles
    for i=1:nSignals
        for j=1:nSignals
            for t=1:nLags
                
                % Statistical sig threshold available for 1 or all lags?
                if nSLags == 1
                    tS = 1;
                else
                    tS = t;
                end
                
                n=min(nBinVect(i),nBinVect(j));
                InormByDist(i,j,t,f)=I(i,j,t,f)/log2(n);
                TnormByDist(i,j,t,f)=T(i,j,t,f)/log2(n);
                Ic(i,j,t,f)=0.5*(1+erf((I(i,j,t,f)-meanShuffI(i,j,tS,f))/(sqrt(2)*sigmaShuffI(i,j,tS,f))));
                Tc(i,j,t,f)=0.5*(1+erf((T(i,j,t,f)-meanShuffT(i,j,tS,f))/(sqrt(2)*sigmaShuffT(i,j,tS,f))));
                SigThreshInormByDist(i,j,f)=SigThreshI(i,j,tS,f)/log2(n);
                SigThreshTnormByDist(i,j,f)=SigThreshT(i,j,tS,f)/log2(n);
                TvsIzero(i,j,t,f)=T(i,j,t,f)/I(i,j,l0i,f);
                SigThreshTvsIzero(i,j,f)=SigThreshT(i,j,tS,f)/I(i,j,l0i,f);
                RelI(i,j,t,f)=I(i,j,t,f)/HYf(i,j,t,f);
                RelT(i,j,t,f)=T(i,j,t,f)/HYf(i,j,t,f);
                HXtNormByDist(i,j,t,f)=HXt(i,j,t,f)/log2(n);
                IvsIzero(i,j,t,f)=I(i,j,t,f)/I(i,j,l0i,f);
                SigThreshIvsIzero(i,j,f)=SigThreshI(i,j,tS,f)/I(i,j,l0i,f);
            end
        end
    end
end

%-----------------------------------------------------------

