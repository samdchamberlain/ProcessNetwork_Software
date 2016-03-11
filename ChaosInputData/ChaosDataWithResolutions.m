%Ben Ruddell - UIUC - 08 October 2007
%A two-way lagged coupled logistic map oscillator with smoothing and time lags both ways

clear all;
close all;
clc;

%SETTINGS
timesteps=10000
lagAB=3
lagBA=3
resAB=10
resBA=10
r=3.99 %this is in the range of chaotic variation 3.57<r<4 for logistic map

%COMPUTATION
for t=1:timesteps
    
    Data(t,1)=rand(1);
    Data(t,2)=rand(1);
    
    if t>lagAB+resAB
        i=t-lagAB;
        i0=i-resAB+1;
        x=mean(Data(i0:i,1)); %map x to the lagged scaled A
        Data(t,2)=r*x*(1-x); %then overwrite B with a logistic map
    end
    
    if t>lagBA+resBA
        i=t-lagBA;
        i0=i-resBA+1;
        x=mean(Data(i0:i,2)); %map x to the lagged scaled B
        Data(t,1)=r*x*(1-x); %then overwrite A with a logistic map        
    end
    
end

%PLOTS
figure(1) %data
plot(Data)
title('The Data')
figure(2) %attractor
title('The Attractor')
scatter(Data(1:timesteps-lagBA,1),Data(lagAB+1:timesteps,2)) 

save data_ChaosWithResolutions Data

%INTERESTING PARAMETER SETS TO TRY...

%1,1,1,1 gives standard-looking logistic map
%1,1,1,10 gives increasing oscillations which then stabilize
%1,1,1,100 gives damped out oscillations which disappear
%1,1,10,10 gives smooth repeating oscillations
%2,2,10,10 expands that attractor a bit
%3,3,10,10 expands that attractor a bit more
%10,10,10,10 gives "slow" smooth chaotic oscillations that synchronize then diverge
%10,10,1,1 gives a choppy chaotic oscillation that synchronizes then diverges on a longer time interval i.e. 10
%100,100,1,1 looks a lot like standard logistic map
%1,100,1,100 gives the smooth chaos again, with scale of about 100
