%Ben Ruddell - UIUC - Noisy Data Generator - 19JAN2007

clear all;
close all;
clc;

disp('Computing an artificial dataset with noise...')

timesteps=10000
sinperiodlength=48
%sineoffset=15
ARcoeffAB1=0.5; %AB1, AB2 must not add up to more than 1....
ARcoeffAB2=0.0;
ARcoeffBA1=0.5; %BA1, BA2 must not add up to more than 1....
ARcoeffBA2=0.0;
ARlagAB1=1;
ARlagAB2=10; %must be more than ARlagAB1
ARlagBA1=7;
ARlagBA2=18; %must be more than ARlagBA1
chaoslagBA=7; %means that the A signal is this many, -1, timesteps behind B
chaoslagAB=1; %means that the B signal is this many timesteps behind A
r=3.99 %this is in the range of chaotic variation 3.57<r<4 for logistic map

NTS = [0 1/100 1/10 1/5 1/2 1/1 2/1 5/1 10/1 100/1 inf] %Noise To Signal Ratios - signal is AR or lag, noise is sin or normal
msize=size(NTS);
nSNRs=msize(2)

nSignals=8

nSin=1
nLagSin=2
nARnoiseA=3
nARnoiseB=4
nNormNoiseA=5
nNormNoiseB=6
nLnoiseA=7
nLnoiseB=8

%generate extra AR noise to allow correlation structure to develop
nARsteps=100000;    %must be greater than "timesteps"
AR=NaN(nARsteps,2);
for t=1:nARsteps
    AR(t,1)=randn(1);
    AR(t,2)=randn(1);
    if t>ARlagBA2
        AR(t,1)=randn(1)+ARcoeffBA1*AR(t-ARlagBA1,2)+ARcoeffBA2*AR(t-ARlagBA2,2);  
    end
    if t>ARlagAB2
        AR(t,2)=randn(1)+ARcoeffAB1*AR(t-ARlagAB1,1)+ARcoeffAB2*AR(t-ARlagAB2,1);    
    end
end
StartStep=nARsteps-timesteps;

%fill in the data
signal=NaN(timesteps,nSignals);
for t=1:timesteps

    %make a sine wave of "timesteps" length that repeats every "period"
    %timesteps and another that's lagged
    signal(t,nSin)=sin(t*2*pi/sinperiodlength);
    %signal(t,nLagSin)=sin((t-sineoffset)*2*pi/sinperiodlength);

    %fill normal noise
    signal(t,nNormNoiseA)=randn(1);
    signal(t,nNormNoiseB)=randn(1);

    %fill AR noise
    t2=t+StartStep-1;
    signal(t,nARnoiseA)=AR(t2,1);
    signal(t,nARnoiseB)=AR(t2,2);
    
    %generate logistic map noise- it's cross-related with the other
    %logistic map across LlagA and LlagB time lags...
    signal(t,nLnoiseA)=rand(1);
    signal(t,nLnoiseB)=rand(1);
    if t>chaoslagBA
        x=signal(t-chaoslagBA,nLnoiseB); %reference the x to the B logistic signal
        signal(t,nLnoiseA)=r*x*(1-x); %then overwrite A with a logistic map
    end
    if t>chaoslagAB
        x=signal(t-chaoslagAB,nLnoiseA); %reference the x to the B logistic signal
        signal(t,nLnoiseB)=r*x*(1-x); %then overwrite A with a logistic map
    end
    
end

% disp('Original Means and Standard Deviations...')
% mean(signal)
% std(signal)

%standardize and normalize all the data to have mean of zero and sigma of 1
means=mean(signal);
sigmas=std(signal);
Usignal=NaN(timesteps,nSignals);
for s=1:nSignals
    for t=1:timesteps
        Usignal(t,s)=(signal(t,s)-means(s))/sigmas(s);
    end
end

% disp('Standardized Means and Standard Deviations... should be zero and one')
% mean(Usignal)
% std(Usignal)

%add portions of the three kinds of noise
ARNoisyData=NaN(timesteps,nSNRs*2);
ARNoisyData2=NaN(timesteps,nSNRs*2);
NormNoisyData=NaN(timesteps,nSNRs*2);
ChaosNoisyData=NaN(timesteps,nSNRs*2);
ChaosNoisyData2=NaN(timesteps,nSNRs*2);
for t=1:timesteps
    index=1;
    for i=1:nSNRs
        
        %AR noise, A and B, for normal and lagged sine
        ARNoisyData(t,index)=Usignal(t,nSin)+Usignal(t,nARnoiseA)*NTS(i);
        ARNoisyData(t,index+1)=Usignal(t,nSin)+Usignal(t,nARnoiseB)*NTS(i);
        ARNoisyData2(t,index)=Usignal(t,nNormNoiseA)+Usignal(t,nARnoiseA)*NTS(i);
        ARNoisyData2(t,index+1)=Usignal(t,nNormNoiseB)+Usignal(t,nARnoiseB)*NTS(i);        

        %normal noise, A and B, for normal and lagged sine
        NormNoisyData(t,index)=Usignal(t,nSin)+Usignal(t,nNormNoiseA)*NTS(i);
        NormNoisyData(t,index+1)=Usignal(t,nSin)+Usignal(t,nNormNoiseB)*NTS(i);

        %Chaotic noisefor normal and lagged sine
        ChaosNoisyData(t,index)=Usignal(t,nSin)+Usignal(t,nLnoiseA)*NTS(i);
        ChaosNoisyData(t,index+1)=Usignal(t,nSin)+Usignal(t,nLnoiseB)*NTS(i);
        ChaosNoisyData2(t,index)=Usignal(t,nNormNoiseA)+Usignal(t,nLnoiseA)*NTS(i);
        ChaosNoisyData2(t,index+1)=Usignal(t,nNormNoiseB)+Usignal(t,nLnoiseB)*NTS(i);        

        %do the noise-only portion as last two indices of each dataset
        if i == nSNRs 

            %AR noise with added corrupting sine
            ARNoisyData(t,index)=Usignal(t,nARnoiseA);
            ARNoisyData(t,index+1)=Usignal(t,nARnoiseB);
            
            %AR noise with added corrupting normal noise
            ARNoisyData2(t,index)=Usignal(t,nARnoiseA);
            ARNoisyData2(t,index+1)=Usignal(t,nARnoiseB);            

            %normal noise with added sine
            NormNoisyData(t,index)=Usignal(t,nNormNoiseA);
            NormNoisyData(t,index+1)=Usignal(t,nNormNoiseB);

            %Chaotic noise with added corrupting sine
            ChaosNoisyData(t,index)=Usignal(t,nLnoiseA);
            ChaosNoisyData(t,index+1)=Usignal(t,nLnoiseB);
            
            %Chaotic noise with added corrupting normal noise
            ChaosNoisyData2(t,index)=Usignal(t,nLnoiseA);
            ChaosNoisyData2(t,index+1)=Usignal(t,nLnoiseB);            

        end

        index=index+2;
        
    end
end

save data_ARNoisy.txt ARNoisyData -ASCII
save data_NormNoisy.txt NormNoisyData -ASCII
save data_ChaosNoisy.txt ChaosNoisyData -ASCII
save data_ARNoisy2.txt ARNoisyData2 -ASCII
save data_ChaosNoisy2.txt ChaosNoisyData2 -ASCII

disp('Plotting the data...')
f=1;
figure(f)
plot(signal)
xlabel('Original Data and Noise')
f=f+1;
figure(f)
plot(Usignal)
xlabel('Standardized Data and Noise')
f=f+1;
figure(f)
plot(ARNoisyData)
xlabel('Proportions of AR noise on Sin Waves')
f=f+1;
figure(f)
plot(NormNoisyData)
xlabel('Proportions of normal noise on Sin Waves')
f=f+1;
figure(f)
plot(ChaosNoisyData)
xlabel('Proportions of logistic noise on Sin Waves')
f=f+1;
figure(f)
plot(ARNoisyData2)
xlabel('Proportions of AR noise with corrupting normal noise')
f=f+1;
figure(f)
plot(ChaosNoisyData2)
xlabel('Proportions of logistic noise with corrupting normal noise')
f=f+1;
figure(f)
hist(signal(:,3),100)
xlabel('histogram of AR noise')
f=f+1;
figure(f)
hist(signal(:,5),100)
xlabel('histogram of normal noise')
f=f+1;
figure(f)
hist(signal(:,7),100)
xlabel('histogram of logistic noise')
