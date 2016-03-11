
function [Abinary,Awtd,AwtdCut,charLagFirstPeak,TcharLagFirstPeak,charLagMaxPeak,TcharLagMaxPeak,TvsIzerocharLagMaxPeak,nSigLags,FirstSigLag,LastSigLag]=AdjMatrices(T,SigThreshT,TvsIzero,lagVect)

% COMPUTES AN ADJACENCY MATRIX A
% COMPUTES THE CHARACTERISTIC LAG WHICH IS THE FIRST SIGNIFICANT LAG
% TAKES THE TRANSFER INFORMATION MATRIX AND THE SIGNIFICANCE THRESHOLDS
% NOTE - This code considers ONLY POSITIVE LAGS in the computation of these matrices

[nSignals,~,nLags,nFiles]=size(T);
nSLags = size(SigThreshT,3);

Abinary=zeros(nSignals,nSignals,nLags,nFiles);
Awtd=NaN(nSignals,nSignals,nLags,nFiles);
AwtdCut=zeros(nSignals,nSignals,nLags,nFiles);
charLagFirstPeak=zeros(nSignals,nSignals,nFiles);
TcharLagFirstPeak=zeros(nSignals,nSignals,nFiles);
charLagMaxPeak=zeros(nSignals,nSignals,nFiles);
TcharLagMaxPeak=zeros(nSignals,nSignals,nFiles);
TvsIzerocharLagMaxPeak=zeros(nSignals,nSignals,nFiles);
nSigLags=zeros(nSignals,nSignals,nFiles);
FirstSigLag=NaN(nSignals,nSignals,nFiles);
LastSigLag=NaN(nSignals,nSignals,nFiles);

% Find index of lag = 0
l0i = find(lagVect == 0,1);

for f=1:nFiles
    for sX=1:nSignals
        for sY=1:nSignals

            FirstPeakFlag=0;
            FirstSigFlag=0;

            Awtd=T;

            
            %check 0 lag
            lag=l0i;
            
            % Statistical sig threshold available for 1 or all lags?
            if nSLags == 1
                lagS = 1;
            else
                lagS = lag;
            end
            
            if T(sX,sY,lag,f) > SigThreshT(sX,sY,lagS,f)
                Abinary(sX,sY,lag,f) = 1;
                AwtdCut(sX,sY,lag,f) = T(sX,sY,lag,f);
                LastSigLag(sX,sY,f)=lag;
                nSigLags(sX,sY,f)=nSigLags(sX,sY,f)+1;
                charLagMaxPeak(sX,sY,f)=lag;
                TcharLagMaxPeak(sX,sY,f)=T(sX,sY,lag,f);
                TvsIzerocharLagMaxPeak(sX,sY,f)=TvsIzero(sX,sY,lag,f);
                FirstSigFlag=1;
                if nLags > 1
                    if T(sX,sY,lag,f) > T(sX,sY,lag+1,f)
                        charLagFirstPeak(sX,sY,f)=lag;
                        TcharLagFirstPeak(sX,sY,f)=T(sX,sY,lag,f);
                        FirstPeakFlag = 1;
                    end
                else
                    charLagFirstPeak(sX,sY,f)=lag;
                    TcharLagFirstPeak(sX,sY,f)=T(sX,sY,lag,f);
                    FirstPeakFlag = 1;
                end
            end

            %check the other lag
            if nLags > 1
                for lag=l0i+1:nLags-1

                    % Statistical sig threshold available for 1 or all lags?
                    if nSLags == 1
                        lagS = 1;
                    else
                        lagS = lag;
                    end
            
                    if T(sX,sY,lag,f) > SigThreshT(sX,sY,lagS,f)
                        Abinary(sX,sY,lag,f) = 1;
                        AwtdCut(sX,sY,lag,f) = T(sX,sY,lag,f);
                        LastSigLag(sX,sY,f)=lag;
                        nSigLags(sX,sY,f)=nSigLags(sX,sY,f)+1;
                        if FirstSigFlag == 0
                            FirstSigLag(sX,sY,f)=lag;
                            FirstSigFlag = 1;
                        end
                        if FirstPeakFlag == 0 && T(sX,sY,lag,f) > T(sX,sY,lag-1,f) && T(sX,sY,lag,f) > T(sX,sY,lag+1,f)
                            charLagFirstPeak(sX,sY,f)=lag;
                            TcharLagFirstPeak(sX,sY,f)=T(sX,sY,lag,f);
                            FirstPeakFlag = 1;
                        end
                        if T(sX,sY,lag) > TcharLagMaxPeak(sX,sY,f)
                            charLagMaxPeak(sX,sY,f)=lag;
                            TcharLagMaxPeak(sX,sY,f)=T(sX,sY,lag,f);
                            TvsIzerocharLagMaxPeak(sX,sY,f)=TvsIzero(sX,sY,lag,f);
                        end
                    end
                end

                %check the last lag
                lag=nLags;
                
                % Statistical sig threshold available for 1 or all lags?
                if nSLags == 1
                    lagS = 1;
                else
                    lagS = lag;
                end
            
                if T(sX,sY,lag,f) > SigThreshT(sX,sY,lagS,f)
                    Abinary(sX,sY,lag,f) = 1;
                    AwtdCut(sX,sY,lag,f) = T(sX,sY,lag,f);
                    LastSigLag(sX,sY,f)=lag;
                    nSigLags(sX,sY,f)=nSigLags(sX,sY,f)+1;
                    if FirstSigFlag == 0
                        FirstSigLag(sX,sY,f)=lag;
                        FirstSigFlag = 1;
                    end
                    if FirstPeakFlag == 0 && T(sX,sY,lag,f) > T(sX,sY,lag-1,f)
                        charLagFirst(sX,sY,f)=lag;
                        TcharLagFirst(sX,sY,f)=T(sX,sY,lag,f);
                        FirstPeakFlag = 1;
                    end
                    if T(sX,sY,lag,f) > TcharLagMaxPeak(sX,sY,f)
                        charLagMaxPeak(sX,sY,f)=lag;
                        TcharLagMaxPeak(sX,sY,f)=T(sX,sY,lag,f);
                        TvsIzerocharLagMaxPeak(sX,sY,f)=TvsIzero(sX,sY,lag,f);
                    end
                end
            end
        end
    end
end