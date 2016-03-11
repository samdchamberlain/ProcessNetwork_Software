function [HXt, HYw, HYf, HXtYw, HXtYf, HYwYf, HXtYwYf, I, T, nCounts] ...
    = ShannonBitsWrapper(classifiedData, lag, nTuples, nBinMat, lagRange, nYw, NoDataCode)

[~,nSignals] = size(classifiedData);
classifiedData(classifiedData == NoDataCode) = NaN;

HXt             = NaN( nSignals, nSignals );
HYw             = NaN( nSignals, nSignals );
HYf             = NaN( nSignals, nSignals );
HXtYw           = NaN( nSignals, nSignals );
HXtYf           = NaN( nSignals, nSignals );
HYwYf           = NaN( nSignals, nSignals );
HXtYwYf         = NaN( nSignals, nSignals );
I               = NaN( nSignals, nSignals );
T               = NaN( nSignals, nSignals );
nCounts         = NaN( nSignals, nSignals );

for sX=1:nSignals
    for sY=1:nSignals
        
        %CONSTRUCT THREE-COLUMN MATRIX WITH COLUMNS TIME-SHIFTED
        XtSTART = max([lagRange(2) nYw])+1-lag;
        YwSTART = max([lagRange(2) nYw])+1-nYw;
        YfSTART = max([lagRange(2) nYw])+1; 

        tupleMat=NaN(nTuples,3);
        tupleMat(:,1)=classifiedData(XtSTART:XtSTART+nTuples-1,sX);        %Leading Node Xt (lag tau earlier than present)
        tupleMat(:,2)=classifiedData(YwSTART:YwSTART+nTuples-1,sY);        %Led Node Yw (present time)
        tupleMat(:,3)=classifiedData(YfSTART:YfSTART+nTuples-1,sY);        %Led Node Yf (one timestep in future)

        %CHECK TO ENSURE TUPLEMAT HAS AT LEAST ONE COMPLETE ROW OF DATA
        if sum(sum(isnan(tupleMat),2) > 0) == nTuples
            logwrite(['Warning: no data in tupleMat, skipping sX = ' num2str(sX) ', sY = ' num2str(sY) ', lag = ' num2str(lag)],1)
            continue
        end
 
        %CALCULATE ENTROPIES FROM TUPLEMAT
        [C, nCounts(sX,sY)] = getCountMat( tupleMat, nBinMat, sX, sY, NaN);
        [HXt(sX,sY),HYw(sX,sY),HYf(sX,sY),HXtYw(sX,sY),HXtYf(sX,sY),HYwYf(sX,sY),HXtYwYf(sX,sY)]=GetShannonBits( C, nCounts(sX,sY) );
        I(sX,sY) = HXt(sX,sY) + HYf(sX,sY) - HXtYf(sX,sY);
        T(sX,sY) = HXtYw(sX,sY) + HYwYf(sX,sY) - HYw(sX,sY) - HXtYwYf(sX,sY);

    end
end