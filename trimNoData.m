function [trimmedData,skipcounter]=trimNoData(rawData,NoDataCode)

[nData,nSignals]=size(rawData);
trimmedData=NaN(nData,nSignals);
NoDataFlag=zeros(nData,1);
skipcounter=0;

for n=1:nData
    for s=1:nSignals
        if rawData(n,s)==NoDataCode
            NoDataFlag(n)=1;
            skipcounter=skipcounter+1;
            break
        end
    end
end

for n2=1:nData
    if NoDataFlag(n2)
        for s=1:nSignals
            trimmedData(n2,s)=NoDataCode;
        end
    else
        for s=1:nSignals
            trimmedData(n2,s)=rawData(n2,s);
        end
    end
end