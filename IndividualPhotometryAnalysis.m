%% Extract photometry data from TDT Data Tank
%UPDATED Jan 2022 for Lerner Lab

%point to folder where data is located
path_to_data = '/Users/jls2314/Documents/MATLAB';

% set up Tank name, variables to extract
tankdir = path_to_data;
tankname = 'dLight-Cohort1';
blockname = 'Photo_140_306-190904-095355';
storenames = {'Dv1A' 'Dv2A' 'Dv3B' 'Dv4B' 'LNRW' 'LNnR' 'PrtR' 'PrtN' 'Sock'}; % name of stores to extract from TDT (usu. 4-letter code)
% extract
for k = 1:numel(storenames)
  storename = storenames{k};
  S{k} = tdt2mat(tankdir, tankname, blockname, storename);
end

%% Extract and massage data

%Get control and signal data structures
Dv1A = S{1}; %Sensor A control signal
Dv2A = S{2}; %Sensor A GCaMP signal
Dv3B = S{3}; %Sensor B control signal
Dv4B = S{4}; %Sensor B GCaMP signal

% Get control and signal data from structures above as a vector (repeat for each channel)
controlA = Dv1A.data;
controlA = reshape(controlA', [],1); % unwrap data from m x 256 array
signalA = Dv2A.data;
signalA = reshape(signalA', [],1); % unwrap data from m x 256 array
controlB = Dv3B.data;
controlB = reshape(controlB', [],1); % unwrap data from m x 256 array
signalB = Dv4B.data;
signalB = reshape(signalB', [],1); % unwrap data from m x 256 array

% Get data timestamps (same for all channels, so just use one)
ts = Dv1A.timestamps;
t_rec_start = ts(1);

ts = ts-ts(1); % convert from Unix time to 'seconds from block start'
ts = bsxfun(@plus, ts(:), (0:Dv1A.npoints-1)*(1./Dv1A.sampling_rate));
ts = reshape(ts',[],1);

%% Get TTL timestamps

% Get timestamps
Control_ts = S{1,5}.timestamps - t_rec_start;
NosepokeRew_ts = S{1,5}.timestamps - t_rec_start;
NosepokeNoRew_ts = S{1,6}.timestamps - t_rec_start;
PortReward_ts = S{1,7}.timestamps - t_rec_start;
PortNoReward_ts = S{1,8}.timestamps - t_rec_start;
Shock_ts = S{1,9}.timestamps - t_rec_start;

%% cut off 1st second (lights turning on)
correctionIndex = ts <1;

controlA(correctionIndex) = [];
signalA(correctionIndex) = [];
ts(correctionIndex) = []; 
controlB(correctionIndex) = [];
signalB(correctionIndex) = [];

%% check data by plotting

figure('Name','SensorA')
hold on;
h1 = plot(ts, controlA(:), 'b');
h2 = plot(ts, signalA(:), 'r');
legend ([h1, h2], 'control A','signal A', 'orientation', 'vertical', 'Location', 'NorthEast');
xlabel('Time (s)');
ylabel('Fluorescence (V from photodetector)');
hold off;


figure('Name','SensorB')
hold on;
h1 = plot(ts, controlB(:), 'b');
h2 = plot(ts, signalB(:), 'r');
legend ([h1, h2], 'control B','signal B', 'orientation', 'vertical', 'Location', 'NorthEast');
xlabel('Time (s)');
ylabel('Fluorescence (V from photodetector)');
hold off;

%% smooth signals
 
controlA = filtfilt(ones(1,100)/100,1, controlA);
signalA = filtfilt(ones(1,100)/100,1, signalA);
controlB = filtfilt(ones(1,100)/100,1, controlB);
signalB = filtfilt(ones(1,100)/100,1, signalB);

%% Fit control to signal and plot to check

[controlFitA] = controlFit (controlA, signalA);
[controlFitB] = controlFit (controlB, signalB);

figure('Name','ControlFit A')
hold on;
h3 = plot(ts, signalA(:), 'r');
h1 = plot (ts, controlFitA(:), 'c');
h2 = plot(ts, controlA(:), 'b');
legend ([h1, h2, h3], 'controlFit A','control A','signal A', 'orientation', 'vertical', 'Location', 'NorthEast');
hold off;

figure('Name','ControlFit B')
hold on;
h3 = plot(ts, signalB(:), 'r');
h1 = plot (ts, controlFitB(:), 'c');
h2 = plot(ts, controlB(:), 'b');
legend ([h1, h2, h3], 'controlFit B','control B','signal B', 'orientation', 'vertical', 'Location', 'NorthEast');
hold off;

%% Get delta F/F using the fitted control signal
 
[normDatA] = deltaFF (signalA, controlFitA);
[normDatB] = deltaFF (signalB, controlFitB);

figure('Name','Normalized Data A')
plot (ts, normDatA, 'k');
figure('Name','Normalized Data B')
plot (ts, normDatB, 'k');

%% Look at where event timestamps are

% plot data on sensor A with timestamps
figure('Name','Data on A with Timestamps')
hold on;
plot (ts, normDatA, 'k');
nControl = numel(Control_ts);
plot(repmat(Control_ts',2,1), repmat([1.2*min(normDatA(:)); 1.2.*max(normDatA(:))], 1, nControl), 'r--');
nNosepokeRew = numel(NosepokeRew_ts);
plot(repmat(NosepokeRew_ts',2,1), repmat([1.2*min(normDatA(:)); 1.2.*max(normDatA(:))], 1, nNosepokeRew), 'g--');
nNosepokeNoRew = numel(NosepokeNoRew_ts);
plot(repmat(NosepokeNoRew_ts',2,1), repmat([1.2*min(normDatA(:)); 1.2.*max(normDatA(:))], 1, nNosepokeNoRew), 'y--');
nPortReward = numel(PortReward_ts);
plot(repmat(PortReward_ts',2,1), repmat([1.2*min(normDatA(:)); 1.2.*max(normDatA(:))], 1, nPortReward), 'o--');
nPortNoReward = numel(PortNoReward_ts);
plot(repmat(PortNoReward_ts',2,1), repmat([1.2*min(normDatA(:)); 1.2.*max(normDatA(:))], 1, nPortNoReward), 'm--');
nShock = numel(Shock_ts);
plot(repmat(Shock_ts',2,1), repmat([1.2*min(normDatA(:)); 1.2.*max(normDatA(:))], 1, nShock), 'o--');


% plot data on sensor B with timestamps
figure('Name','Data on B with Timestamps')
hold on;
plot (ts, normDatB, 'k');
nControl = numel(Control_ts);
plot(repmat(Control_ts',2,1), repmat([1.2*min(normDatB(:)); 1.2.*max(normDatB(:))], 1, nControl), 'r--');
nNosepokeRew = numel(NosepokeRew_ts);
plot(repmat(NosepokeRew_ts',2,1), repmat([1.2*min(normDatB(:)); 1.2.*max(normDatB(:))], 1, nNosepokeRew), 'g--');
nNosepokeNoRew = numel(NosepokeNoRew_ts);
plot(repmat(NosepokeNoRew_ts',2,1), repmat([1.2*min(normDatB(:)); 1.2.*max(normDatB(:))], 1, nNosepokeNoRew), 'y--');
nPortReward = numel(PortReward_ts);
plot(repmat(PortReward_ts',2,1), repmat([1.2*min(normDatB(:)); 1.2.*max(normDatB(:))], 1, nPortReward), 'b--');
nPortNoReward = numel(PortNoReward_ts);
plot(repmat(PortNoReward_ts',2,1), repmat([1.2*min(normDatB(:)); 1.2.*max(normDatB(:))], 1, nPortNoReward), 'm--');
nShock = numel(Shock_ts);
plot(repmat(Shock_ts',2,1), repmat([1.2*min(normDatB(:)); 1.2.*max(normDatB(:))], 1, nShock), 'o--');

%% Get z-scores

zDatA = (normDatA - nanmean(normDatA))/ nanstd(normDatA);
zDatB = (normDatB - nanmean(normDatB))/ nanstd(normDatB);

figure('Name','A')
plot (ts, zDatA, 'k');
ylabel('z-score');
figure('Name','B')
plot (ts, zDatB, 'k');
ylabel('z-score');

%% make PSTH arrays based on reward/press timestamps, then average rows into a mean vector for plotting
 
samplingRate = Dv1A.sampling_rate; % TDT returns 'true' sample rate. Use any channel.

%Define window for PSTH, seconds before and after to look
nSecPrev = 10;
nSecPost = 20; 
% convert seconds to TDT timestamps
nTsPrev = round (nSecPrev * samplingRate);
nTsPost = round (nSecPost * samplingRate);


%% further timestamp processing - delete timestamps beyond the first from bursts of pressing for 
%Right Non rewarded nosepokes
nNosepokeNoRew_ts = size(NosepokeNoRew_ts,1);
NosepokeNoRewextra_ts = NaN(nNosepokeNoRew_ts,1);
for iNosepokeNoRewextra = 2:nNosepokeNoRew_ts-1;
    thisTime = NosepokeNoRew_ts(iNosepokeNoRewextra);
    prevTime = NosepokeNoRew_ts(iNosepokeNoRewextra-1);
    if (thisTime - prevTime) < 2
        NosepokeNoRewextra_ts (iNosepokeNoRewextra) = thisTime;
    end
end
% Make NoRewd_ts and delete "extras"
NoBurstNosepokeNoRew_ts = setdiff (NosepokeNoRew_ts, NosepokeNoRewextra_ts);
%% further timestamp processing - delete timestamps beyond the first from bursts of pressing or mag entry
 % Make a list of Non-Rewarded Port Entries that are larger than 1 sec
 % apart
nPortNoReward_ts = size(PortNoReward_ts,1);
Magextra_ts = NaN(nPortNoReward_ts,1);
for iMagNo = 2:nPortNoReward_ts-1;
    thisTime = PortNoReward_ts(iMagNo);
    prevTime = PortNoReward_ts(iMagNo-1);
    if (thisTime - prevTime) < 2
        Magextra_ts (iMagNo) = thisTime;
    end
end
% Make NoRewd_ts and delete "extras"
NoBurstPortNoReward_ts = setdiff (PortNoReward_ts, Magextra_ts);

%%
%if plotting Inactive Nosepoke
nControl_ts = size(Control_ts,1);
%A
ControlPsthArrayA = NaN(nControl_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iControl = 1:nControl_ts
    thisTime = Control_ts(iControl)-1;
    thisIndex = round((thisTime*samplingRate));
    ControlPsthArrayA(iControl,:) = processPhotDataRow_normDat(zDatA, thisIndex, nTsPrev, nTsPost);
end
ControlErrA = (nanstd(ControlPsthArrayA))/sqrt(size(ControlPsthArrayA,1));
ControlPsthA = nanmean(ControlPsthArrayA);
%B
ControlPsthArrayB = NaN(nControl_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iControl = 1:nControl_ts
    thisTime = Control_ts(iControl)-1;
    thisIndex = round((thisTime*samplingRate));
    ControlPsthArrayB(iControl,:) = processPhotDataRow_normDat(zDatB, thisIndex, nTsPrev, nTsPost);
end
ControlErrB = (nanstd(ControlPsthArrayB))/sqrt(size(ControlPsthArrayB,1));
ControlPsthB = nanmean(ControlPsthArrayB,1);
%% if plotting Rewarded Nosepoke
nNosepokeRew_ts = size(NosepokeRew_ts,1);
%A
NosepokeRewPsthArrayA = NaN(nNosepokeRew_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iNosepokeRew = 1:nNosepokeRew_ts
    thisTime = NosepokeRew_ts(iNosepokeRew)-1;
    thisIndex = round((thisTime*samplingRate));
    NosepokeRewPsthArrayA(iNosepokeRew,:) = processPhotDataRow_normDat(zDatA, thisIndex, nTsPrev, nTsPost);
end
NosepokeRewErrA = (nanstd(NosepokeRewPsthArrayA))/sqrt(size(NosepokeRewPsthArrayA,1));
NosepokeRewPsthA = nanmean(NosepokeRewPsthArrayA);
%B
NosepokeRewPsthArrayB = NaN(nNosepokeRew_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iNosepokeRew = 1:nNosepokeRew_ts
    thisTime = NosepokeRew_ts(iNosepokeRew)-1;
    thisIndex = round((thisTime*samplingRate));
    NosepokeRewPsthArrayB(iNosepokeRew,:) = processPhotDataRow_normDat(zDatB, thisIndex, nTsPrev, nTsPost);
end
NosepokeRewErrB = (nanstd(NosepokeRewPsthArrayB))/sqrt(size(NosepokeRewPsthArrayB,1));
NosepokeRewPsthB = nanmean(NosepokeRewPsthArrayB);


%% if plotting Unrewarded Nosepoke
% if cutting off time, change to:
%nNoBurstNosepokeNoRew_ts= number of timestamps remaining-1
nNoBurstNosepokeNoRew_ts = size(NoBurstNosepokeNoRew_ts,1);
%A
NosepokeNoRewPsthArrayA = NaN(nNoBurstNosepokeNoRew_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iNoBurstNosepokeNoRew = 1:nNoBurstNosepokeNoRew_ts
    thisTime = NoBurstNosepokeNoRew_ts(iNoBurstNosepokeNoRew)-1;
    thisIndex = round((thisTime*samplingRate));
    NosepokeNoRewPsthArrayA(iNoBurstNosepokeNoRew,:) = processPhotDataRow_normDat(zDatA, thisIndex, nTsPrev, nTsPost);
end
NosepokeNoRewErrA = (nanstd(NosepokeNoRewPsthArrayA))/sqrt(size(NosepokeNoRewPsthArrayA,1));
NosepokeNoRewPsthA = nanmean(NosepokeNoRewPsthArrayA);
%B
NosepokeNoRewPsthArrayB = NaN(nNoBurstNosepokeNoRew_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iNoBurstNosepokeNoRew = 1:nNoBurstNosepokeNoRew_ts
    thisTime = NoBurstNosepokeNoRew_ts(iNoBurstNosepokeNoRew)-1;
    thisIndex = round((thisTime*samplingRate));
    NosepokeNoRewPsthArrayB(iNoBurstNosepokeNoRew,:) = processPhotDataRow_normDat(zDatB, thisIndex, nTsPrev, nTsPost);
end
NosepokeNoRewErrB = (nanstd(NosepokeNoRewPsthArrayB))/sqrt(size(NosepokeNoRewPsthArrayB,1));
NosepokeNoRewPsthB = nanmean(NosepokeNoRewPsthArrayB);

%% if plotting PortReward
nPortReward_ts = size(PortReward_ts,1);
%A
PortRewardPsthArrayA = NaN(nPortReward_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iPortReward = 1:nPortReward_ts
    thisTime = PortReward_ts(iPortReward)-1;
    thisIndex = round((thisTime*samplingRate));
    PortRewardPsthArrayA(iPortReward,:) = processPhotDataRow_normDat(zDatA, thisIndex, nTsPrev, nTsPost);
end
PortRewardErrA = (nanstd(PortRewardPsthArrayA))/sqrt(size(PortRewardPsthArrayA,1));
PortRewardPsthA = nanmean(PortRewardPsthArrayA);
%B
PortRewardPsthArrayB = NaN(nPortReward_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iPortReward = 1:nPortReward_ts
    thisTime = PortReward_ts(iPortReward)-1;
    thisIndex = round((thisTime*samplingRate));
    PortRewardPsthArrayB(iPortReward,:) = processPhotDataRow_normDat(zDatB, thisIndex, nTsPrev, nTsPost);
end
PortRewardErrB = (nanstd(PortRewardPsthArrayB))/sqrt(size(PortRewardPsthArrayB,1));
PortRewardPsthB = nanmean(PortRewardPsthArrayB);
% 
%% if plotting PortNoReward
nNoBurstPortNoReward_ts = size(NoBurstPortNoReward_ts,1);
%A
PortNoRewardPsthArrayA = NaN(nNoBurstPortNoReward_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iNoBurstPortNoReward = 1:nNoBurstPortNoReward_ts
    thisTime = NoBurstPortNoReward_ts(iNoBurstPortNoReward)-1;
    thisIndex = round((thisTime*samplingRate));
    PortNoRewardPsthArrayA(iNoBurstPortNoReward,:) = processPhotDataRow_normDat(zDatA, thisIndex, nTsPrev, nTsPost);
end
PortNoRewardErrA = (nanstd(PortNoRewardPsthArrayA))/sqrt(size(PortNoRewardPsthArrayA,1));
PortNoRewardPsthA = nanmean(PortNoRewardPsthArrayA);
%B
PortNoRewardPsthArrayB = NaN(nNoBurstPortNoReward_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iNoBurstPortNoReward = 1:nNoBurstPortNoReward_ts
    thisTime = NoBurstPortNoReward_ts(iNoBurstPortNoReward)-1;
    thisIndex = round((thisTime*samplingRate));
    PortNoRewardPsthArrayB(iNoBurstPortNoReward,:) = processPhotDataRow_normDat(zDatB, thisIndex, nTsPrev, nTsPost);
end
PortNoRewardErrB = (nanstd(PortNoRewardPsthArrayB))/sqrt(size(PortNoRewardPsthArrayB,1));
PortNoRewardPsthB = nanmean(PortNoRewardPsthArrayB);

%% For plotting Shocks
nShock_ts = size(Shock_ts,1);
% A
ShockPsthArrayA = NaN(nShock_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iShock = 1:nShock_ts
    thisTime = Shock_ts(iShock)-1;
    thisIndex = round((thisTime*samplingRate));
    ShockPsthArrayA(iShock,:) = processPhotDataRow_normDat(zDatA, thisIndex, nTsPrev, nTsPost);
end
ShockErrA = (nanstd(ShockPsthArrayA))/sqrt(size(ShockPsthArrayA,1));
ShockPsthA = nanmean(ShockPsthArrayA);
% B
ShockPsthArrayB = NaN(nShock_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iShock = 1:nShock_ts
    thisTime = Shock_ts(iShock)-1;
    thisIndex = round((thisTime*samplingRate));
    ShockPsthArrayB(iShock,:) = processPhotDataRow_normDat(zDatB, thisIndex, nTsPrev, nTsPost);
end
ShockErrB = (nanstd(ShockPsthArrayB))/sqrt(size(ShockPsthArrayB,1));
ShockPsthB = nanmean(ShockPsthArrayB);
%%
%%Baseline Adjustment Section 
%Calculating Baseline Average for Signals A and B immediately prior to Rewarded Nosepoke
 for iNosepokeRew = 1:nNosepokeRew
    BaselineNPA(iNosepokeRew,:) = mean(NosepokeRewPsthArrayA(iNosepokeRew,9157:10174));
end
BaselineNPB = NaN(nNosepokeRew:1);
for iNosepokeRew = 1:nNosepokeRew
    BaselineNPB(iNosepokeRew,:) = mean(NosepokeRewPsthArrayB(iNosepokeRew,9157:10174));
end
% Calculating Baseline Average for Signals A and B immediately prior to Unrewarded Nosepoke
BaselineNoNPA = NaN(nNoBurstNosepokeNoRew_ts:1);
for iNoBurstNosepokeNoRew = 1:nNoBurstNosepokeNoRew_ts
    BaselineNoNPA(iNoBurstNosepokeNoRew,:) = mean(NosepokeNoRewPsthArrayA(iNoBurstNosepokeNoRew,9157:10174));
end
BaselineNoNPB = NaN(nNoBurstNosepokeNoRew_ts:1);
for iNoBurstNosepokeNoRew = 1:nNoBurstNosepokeNoRew_ts
    BaselineNoNPB(iNoBurstNosepokeNoRew,:) = mean(NosepokeNoRewPsthArrayB(iNoBurstNosepokeNoRew,9157:10174));
end
% Calculating Baseline Average for Signals A and B immediately prior to Rewarded Port Entry
BaselineRPA = NaN(nPortReward:1);
for iPortReward = 1:nPortReward
    BaselineRPA(iPortReward,:) = mean(PortRewardPsthArrayA(iPortReward,9157:10174));
end
BaselineRPB = NaN(nPortReward:1);
for iPortReward = 1:nPortReward
    BaselineRPB(iPortReward,:) = mean(PortRewardPsthArrayB(iPortReward,9157:10174));
end
% Calculating Baseline Average for Signals A and B immediately prior to Unewarded Port Entry
BaselineNoRPA = NaN(nNoBurstPortNoReward_ts:1);
for iNoBurstPortNoReward = 1:nNoBurstPortNoReward_ts
    BaselineNoRPA(iNoBurstPortNoReward,:) = mean(PortNoRewardPsthArrayA(iNoBurstPortNoReward,9157:10174));
end
BaselineNoRPB = NaN(nNoBurstPortNoReward_ts:1);
for iNoBurstPortNoReward = 1:nNoBurstPortNoReward_ts
    BaselineNoRPB(iNoBurstPortNoReward,:) = mean(PortNoRewardPsthArrayB(iNoBurstPortNoReward,9157:10174));
end

%Calculating Baseline Average for Signals A and B immediately prior to
%Shock
BaselineShockA = NaN(nShock_ts:1);
for iShock = 1:nShock_ts
    BaselineShockA(iShock,:) = mean(ShockPsthArrayA(iShock,9157:10174));
end
BaselineShockB = NaN(nShock_ts:1);
for iShock = 1:nShock_ts
    BaselineShockB(iShock,:) = mean(ShockPsthArrayB(iShock,9157:10174));
end

%% Calculating Adjusted Values for Rewarded Nosepokes
NPAadjusted = NaN(nNosepokeRew,nTsPrev+nTsPost+1);
for dNosepokeRew = 1:nNosepokeRew
    NPAadjusted(dNosepokeRew,:) = ((NosepokeRewPsthArrayA(dNosepokeRew,1:30519))-(BaselineNPA(dNosepokeRew,:)));
end
NosepokeADJRewErrA = (nanstd(NPAadjusted))/sqrt(size(NPAadjusted,1));
NosepokeADJRewPsthA = nanmean(NPAadjusted);

NPBadjusted = NaN(nNosepokeRew,nTsPrev+nTsPost+1);
for dNosepokeRew = 1:nNosepokeRew
    NPBadjusted(dNosepokeRew,:) = ((NosepokeRewPsthArrayB(dNosepokeRew,1:30519))-(BaselineNPB(dNosepokeRew,:)));
end
NosepokeADJRewErrB = (nanstd(NPBadjusted))/sqrt(size(NPBadjusted,1));
NosepokeADJRewPsthB = nanmean(NPBadjusted);

%Calculating Adjusted Values for Unrewarded Nosepoke
NoNPAadjusted = NaN(nNoBurstNosepokeNoRew_ts,nTsPrev+nTsPost+1);
for dNosepokeNoRew = 1:nNoBurstNosepokeNoRew_ts
    NoNPAadjusted(dNosepokeNoRew,:) = ((NosepokeNoRewPsthArrayA(dNosepokeNoRew,1:30519))-(BaselineNoNPA(dNosepokeNoRew,:)));
end
NoNosepokeADJRewErrA = (nanstd(NoNPAadjusted))/sqrt(size(NoNPAadjusted,1));
NoNosepokeADJRewPsthA = nanmean(NoNPAadjusted);

NoNPBadjusted = NaN(nNoBurstNosepokeNoRew_ts,nTsPrev+nTsPost+1);
for dNosepokeNoRew = 1:nNoBurstNosepokeNoRew_ts
    NoNPBadjusted(dNosepokeNoRew,:) = ((NosepokeNoRewPsthArrayB(dNosepokeNoRew,1:30519))-(BaselineNoNPB(dNosepokeNoRew,:)));
end
NoNosepokeADJRewErrB = (nanstd(NoNPBadjusted))/sqrt(size(NoNPBadjusted,1));
NoNosepokeADJRewPsthB = nanmean(NoNPBadjusted);

% Calculating Adjusted Values for Rewarded Port Entries
RPAadjusted = NaN(nPortReward,nTsPrev+nTsPost+1);
for dPortReward = 1:nPortReward
    RPAadjusted(dPortReward,:) = ((PortRewardPsthArrayA(dPortReward,1:30519))-(BaselineRPA(dPortReward,:)));
end
PortRewardADJRewErrA = (nanstd(RPAadjusted))/sqrt(size(RPAadjusted,1));
PortRewardADJRewPsthA = nanmean(RPAadjusted);

RPBadjusted = NaN(nPortReward,nTsPrev+nTsPost+1);
for dPortReward = 1:nPortReward
    RPBadjusted(dPortReward,:) = ((PortRewardPsthArrayB(dPortReward,1:30519))-(BaselineRPB(dPortReward,:)));
end
PortRewardADJRewErrB = (nanstd(RPBadjusted))/sqrt(size(RPBadjusted,1));
PortRewardADJRewPsthB = nanmean(RPBadjusted);

% Calculating Adjusted Values for Unrewarded Port Entries
NoRPAadjusted = NaN(nNoBurstPortNoReward_ts,nTsPrev+nTsPost+1);
for dPortNoReward = 1:nNoBurstPortNoReward_ts
    NoRPAadjusted(dPortNoReward,:) = ((PortNoRewardPsthArrayA(dPortNoReward,1:30519))-(BaselineNoRPA(dPortNoReward,:)));
end
NoPortRewardADJRewErrA = (nanstd(NoRPAadjusted))/sqrt(size(NoRPAadjusted,1));
NoPortRewardADJRewPsthA = nanmean(NoRPAadjusted);

NoRPBadjusted = NaN(nNoBurstPortNoReward_ts,nTsPrev+nTsPost+1);
for dPortNoReward = 1:nNoBurstPortNoReward_ts
    NoRPBadjusted(dPortNoReward,:) = ((PortNoRewardPsthArrayB(dPortNoReward,1:30519))-(BaselineNoRPB(dPortNoReward,:)));
end
NoPortRewardADJRewErrB = (nanstd(NoRPBadjusted))/sqrt(size(NoRPBadjusted,1));
NoPortRewardADJRewPsthB = nanmean(NoRPBadjusted);

% Calculating Adjusted Values for Shock
ShockAadjusted = NaN(nShock_ts,nTsPrev+nTsPost+1);
for dShock = 1:nShock_ts
    ShockAadjusted(dShock,:) = ((ShockPsthArrayA(dShock,1:30519))-(BaselineShockA(dShock,:)));
end
ShockADJErrA = (nanstd(ShockAadjusted))/sqrt(size(ShockAadjusted,1));
ShockADJPsthA = nanmean(ShockAadjusted);

ShockBadjusted = NaN(nShock_ts,nTsPrev+nTsPost+1);
for dShock = 1:nShock_ts
    ShockBadjusted(dShock,:) = ((ShockPsthArrayB(dShock,1:30519))-(BaselineShockB(dShock,:)));
end
ShockADJErrB = (nanstd(ShockBadjusted))/sqrt(size(ShockBadjusted,1));
ShockADJPsthB = nanmean(ShockBadjusted);


%% Make Plots and find Peaks/AUC
%Rewarded Nosepoke
totalTs = nTsPrev + nTsPost;
increment = (nSecPrev + nSecPost) / totalTs;
timeAxis = (-1 * nSecPrev) : increment : nSecPost;
NosepokeADJRewErrPosA = NosepokeADJRewPsthA + NosepokeADJRewErrA;
NosepokeADJRewErrNegA = NosepokeADJRewPsthA - NosepokeADJRewErrA;
NosepokeADJRewErrPosB = NosepokeADJRewPsthB + NosepokeADJRewErrB;
NosepokeADJRewErrNegB = NosepokeADJRewPsthB - NosepokeADJRewErrB;
figure('Name', 'Rewarded Nosepoke', 'Color',[1 1 1])
hold on
xlim([-1,10]);
set(gca, 'Fontsize', 14);
h1 = plot(timeAxis,NosepokeADJRewPsthA, 'b-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[(NosepokeADJRewErrPosA), fliplr(NosepokeADJRewErrNegA)], [0,0,1],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis,NosepokeADJRewPsthB, 'r-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[(NosepokeADJRewErrPosB), fliplr(NosepokeADJRewErrNegB)], [1,0,0],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontSize', 12);
ylabel('Z-score', 'FontSize', 12);
legend ([h1, h2], 'DMS','DLS', 'orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;
% DMSpeak = max (NosepokeADJRewPsthA)
DMSpeakSeparated = max(NPAadjusted,[],2);   % finds maximums for each individual recording of A
% DLSpeak = max (NosepokeADJRewPsthB)
DLSpeakSeparated = max(NPBadjusted,[],2);   % finds maximums for each individual recording of B

DMSvalley = min (NosepokeRewPsthA);
DMSvalleySeparated = min(NPAadjusted,[],2);   % finds the minimum for each individual recording of A

maxNosepokeRewPsthA = max((NosepokeADJRewPsthA(1,10173:12207)))

maxNosepokeRewPsthB = max((NosepokeADJRewPsthB(1,10173:12207)))

% finding area under curve
x = [0:1:10173];
y = NosepokeADJRewPsthB(1,10173:20346);
DLSareaUC = trapz(x,y)

% finding area under curve for DMS
x = [0:1:10173];
y = NosepokeADJRewPsthA(1,10173:20346);
DMSareaUC = trapz(x,y)


%% make Shock plot
totalTs = nTsPrev + nTsPost;
increment = (nSecPrev + nSecPost) / totalTs;
timeAxis = (-1 * nSecPrev) : increment : nSecPost;
ShockADJErrPosA = ShockADJPsthA + ShockADJErrA;
ShockADJErrNegA = ShockADJPsthA - ShockADJErrA;
ShockADJErrPosB = ShockADJPsthB + ShockADJErrB;
ShockADJErrNegB = ShockADJPsthB - ShockADJErrB;
figure('Name', 'Shock', 'Color',[1 1 1])
hold on
xlim([-1,10]);
set(gca, 'Fontsize', 14);
h1 = plot(timeAxis,ShockADJPsthA, 'b-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[(ShockADJErrPosA), fliplr(ShockADJErrNegA)], [0,0,1],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis,ShockADJPsthB, 'r-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[(ShockADJErrPosB), fliplr(ShockADJErrNegB)], [1,0,0],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontSize', 12);
ylabel('Z-score', 'FontSize', 12);
legend ([h1, h2], 'DMS','DLS', 'orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;
% DMSpeak = max (NosepokeADJRewPsthA)
DMSpeakSeparated = max(ShockAadjusted,[],2);   % finds maximums for each individual recording of A
% DLSpeak = max (NosepokeADJRewPsthB)
DLSpeakSeparated = max(ShockBadjusted,[],2);   % finds maximums for each individual recording of B

DMSvalley = min (ShockPsthA)
DMSvalleySeparated = min(ShockAadjusted,[],2);   % finds the minimum for each individual recording of A
DLSvalley = min (ShockPsthB)
maxShockPsthA = max((ShockADJPsthA(1,10173:12207)))

maxShockPsthB = max((ShockADJPsthB(1,10173:12207)))

% finding area under curve
x = [0:1:10173];
y = ShockADJPsthB(1,10173:20346);
DLSareaUC = trapz(x,y)

% finding area under curve for DMS
x = [0:1:10173];
y = ShockADJPsthA(1,10173:20346);
DMSareaUC = trapz(x,y)
%% make  Nosepoke Unrewarded plot
totalTs = nTsPrev + nTsPost;
increment = (nSecPrev + nSecPost) / totalTs;
timeAxis = (-1 * nSecPrev) : increment : nSecPost;
NoNosepokeADJRewErrPosA = NoNosepokeADJRewPsthA + NoNosepokeADJRewErrA;
NoNosepokeADJRewErrNegA = NoNosepokeADJRewPsthA - NoNosepokeADJRewErrA;
NoNosepokeADJRewErrPosB = NoNosepokeADJRewPsthB + NoNosepokeADJRewErrB;
NoNosepokeADJRewErrNegB = NoNosepokeADJRewPsthB - NoNosepokeADJRewErrB;
figure('Name', 'Unrewarded Nosepoke', 'Color',[1 1 1])
hold on
xlim([-1,10]);
set(gca, 'Fontsize', 14);
h1 = plot(timeAxis,NoNosepokeADJRewPsthA, 'b-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[(NoNosepokeADJRewErrPosA), fliplr(NoNosepokeADJRewErrNegA)], [0,0,1],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis,NoNosepokeADJRewPsthB, 'r-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[(NoNosepokeADJRewErrPosB), fliplr(NoNosepokeADJRewErrNegB)], [1,0,0],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontSize', 12);
ylabel('Z-score', 'FontSize', 12);
legend ([h1, h2], 'DMS','DLS', 'orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;

maxNosepokeNoRewPsthA = max(abs(NoNosepokeADJRewPsthA(1,10173:12207)))

maxNosepokeNoRewPsthB = max(abs(NoNosepokeADJRewPsthB(1,10173:12207)));

%% make  Reward Port Entry plot
totalTs = nTsPrev + nTsPost;
increment = (nSecPrev + nSecPost) / totalTs;
timeAxis = (-1 * nSecPrev) : increment : nSecPost;
PortRewardADJRewErrPosA = PortRewardADJRewPsthA + PortRewardADJRewErrA;
PortRewardADJRewErrNegA = PortRewardADJRewPsthA - PortRewardADJRewErrA;
PortRewardADJRewErrPosB = PortRewardADJRewPsthB + PortRewardADJRewErrB;
PortRewardADJRewErrNegB = PortRewardADJRewPsthB - PortRewardADJRewErrB;
figure('Name', 'Rewarded Port Entry', 'Color',[1 1 1])
hold on
xlim([-1,10]);
set(gca, 'Fontsize', 14);
h1 = plot(timeAxis,PortRewardADJRewPsthA, 'b-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[(PortRewardADJRewErrPosA), fliplr(PortRewardADJRewErrNegA)], [0,0,1],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis,PortRewardADJRewPsthB, 'r-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[(PortRewardADJRewErrPosB), fliplr(PortRewardADJRewErrNegB)], [1,0,0],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontSize', 12);
ylabel('Z-score', 'FontSize', 12);
legend ([h1, h2], 'DMS','DLS', 'orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;

%% make Unreward Port Entry plot
totalTs = nTsPrev + nTsPost;
increment = (nSecPrev + nSecPost) / totalTs;
timeAxis = (-1 * nSecPrev) : increment : nSecPost;
NoPortRewardADJRewErrPosA = NoPortRewardADJRewPsthA + NoPortRewardADJRewErrA;
NoPortRewardADJRewErrNegA = NoPortRewardADJRewPsthA - NoPortRewardADJRewErrA;
NoPortRewardADJRewErrPosB = NoPortRewardADJRewPsthB + NoPortRewardADJRewErrB;
NoPortRewardADJRewErrNegB = NoPortRewardADJRewPsthB - NoPortRewardADJRewErrB;
figure('Name', 'Unrewarded Port Entry', 'Color',[1 1 1])
hold on
xlim([-1,10]);
set(gca, 'Fontsize', 14);
h1 = plot(timeAxis,NoPortRewardADJRewPsthA, 'b-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[(NoPortRewardADJRewErrPosA), fliplr(NoPortRewardADJRewErrNegA)], [0,0,1],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis,NoPortRewardADJRewPsthB, 'r-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[(NoPortRewardADJRewErrPosB), fliplr(NoPortRewardADJRewErrNegB)], [1,0,0],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontSize', 12);
ylabel('Z-score', 'FontSize', 12);
legend ([h1, h2], 'DMS','DLS', 'orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;


%%
%%Behavioral Analysis 
Combine Nosepokes for calculating average NPs before a port entry
CombineNP_ts = [Control_ts; NosepokeRew_ts; NosepokeNoRew_ts];
TotalNP_ts = sort(CombineNP_ts);
CombinePE_ts = [PortReward_ts; PortNoReward_ts];
TotalPE_ts = sort(CombinePE_ts);
    
%% Putting # of NPs before a port entry into an array
nPortReward_ts = size(PortReward_ts,1);
NPsbeforeRPE = NaN(nPortReward_ts, 1);
for iPortReward = 1:nPortReward_ts
    if iPortReward >= 2
        A = sum(TotalNP_ts > PortReward_ts(iPortReward-1,1) & TotalNP_ts < PortReward_ts(iPortReward,1));
        NPsbeforeRPE (iPortReward, 1) = A;
    else
        A = sum(TotalNP_ts < PortReward_ts(iPortReward,1)); 
        NPsbeforeRPE (iPortReward, 1) = A;
    end
end   
AVGNosepokeb4RPE = mean(NPsbeforeRPE);

%% Putting # of NPs before an unrewarded port entry into an array
nPortNoReward_ts = size(PortNoReward_ts,1);
NPsbeforeNPE = NaN(nPortNoReward_ts, 1);
for iPortNoReward = 1:nPortNoReward_ts
    if iPortNoReward >= 2
        B = sum(TotalNP_ts > PortNoReward_ts(iPortNoReward-1,1) & TotalNP_ts < PortNoReward_ts(iPortNoReward,1));
        NPsbeforeNPE (iPortNoReward, 1) = B;
    else
        B = sum(TotalNP_ts < PortNoReward_ts(iPortNoReward,1)); 
        NPsbeforeNPE (iPortNoReward, 1) = B;
    end
end   
AVGNosepokeb4NPE = mean(NPsbeforeNPE);

%% Average time between Port Entries 
AllViablePE=[NoBurstPortNoReward_ts; PortReward_ts];
nAllViablePE = size(AllViablePE,1);
LatencybeforePE = NaN(nAllViablePE,1);
for iAllViablePE = 1:nAllViablePE
    if iAllViablePE >= 2
        Latency = AllViablePE(iAllViablePE,1)- AllViablePE(iAllViablePE-1,1);
        LatencybeforeRPE (iAllViablePE, 1) = Latency;
    else 
        Latency = AllViablePE(iAllViablePE, 1);
        LatencybeforeRPE (iAllViablePE, 1) = Latency;
    end     
end
AVGTimebnPE = mean(LatencybeforeRPE);


%% Average time from NP to Port Reward
nPortReward_ts = size(PortReward_ts,1);
TimeNPtoRPE = NaN(nPortReward_ts,1);
for iPortReward = 1:nPortReward_ts
        AllValues=find(TotalNP_ts < PortReward_ts(iPortReward,1));
        Valueslessthan = TotalNP_ts(AllValues);
        PreviousNP=max(Valueslessthan);
        ValuesLessThanPE (iPortReward,1) = PreviousNP; 
        NPtoRPELatency=PortReward_ts(iPortReward,1) - PreviousNP;
        TimeNPtoRPE(iPortReward, 1) = NPtoRPELatency;
end
AvgNPtoRewardPort = mean(TimeNPtoRPE);

%%
GreatestNP = max(TotalNP_ts);
NoBurstNPEViable1 = find(NoBurstNPE_ts > TotalNP_ts(1,1));
NoBurstNPEViable1_ts = NoBurstNPE_ts (NoBurstNPEViable1);
NoBurstNPEViable = find(NoBurstNPEViable1_ts < GreatestNP);
NoBurstNPEViable_ts = NoBurstNPEViable1_ts (NoBurstNPEViable);
%%
%Create list of Port Entries that directly follow a nosepoke
nNoBurstNPEViable_ts = size(NoBurstNPEViable_ts,1);
IndividualNonRewarded = NaN(nNoBurstNPEViable_ts,1);
for iNoBurstNPEViable = 1:nNoBurstNPEViable_ts
LowerValues=find(TotalNP_ts < NoBurstNPEViable_ts(iNoBurstNPEViable,1));
  if LowerValues>0   
      LowerNPValuesNonRewarded = TotalNP_ts(LowerValues);
      LowestPreviousNP=max(LowerNPValuesNonRewarded);
  end    
UpperValues=find(TotalNP_ts > NoBurstNPEViable_ts(iNoBurstNPEViable,1));
    if UpperValues>0
        UpperNPValuesNonRewarded = TotalNP_ts(UpperValues);
        HighestPreviousNP=min(UpperNPValuesNonRewarded);   
    end    
RangeofPEltoh = find(NoBurstNPEViable_ts > LowestPreviousNP); 
RangeofPElow_ts = NoBurstNPEViable_ts(RangeofPEltoh);
RangeofPEhtol = find(NoBurstNPEViable_ts < HighestPreviousNP);
RangeofPEhigh_ts = NoBurstNPEViable_ts(RangeofPEhtol);
ExtraPortEntries_ts = intersect(RangeofPEhigh_ts, RangeofPElow_ts);
IndividualPoke = min(ExtraPortEntries_ts);
IndividualNonRewarded(iNoBurstNPEViable, 1) = IndividualPoke;        
end
% Make NoRewd_ts and delete "extras"
FinalNonReward_ts = unique(IndividualNonRewarded);

%% Average time between Port Entries - Non Rewarded
nFinalNonReward_ts = size(FinalNonReward_ts,1);
TimeNPtoNPE = NaN(nFinalNonReward_ts,1);
for iFinalNonReward = 1:nFinalNonReward_ts
        AllNONValues=find(TotalNP_ts < FinalNonReward_ts(iFinalNonReward,1));
        NONValueslessthan = TotalNP_ts(AllNONValues);
        NONPreviousNP=max(NONValueslessthan);
        NPtoNPELatency=FinalNonReward_ts(iFinalNonReward,1) - NONPreviousNP;
        TimeNPtoNPE(iFinalNonReward, 1) = NPtoNPELatency;
end
AvgNPtoNoPort = mean(TimeNPtoNPE);

%%
% if plotting RightNosepokeRew
nNosepokeRew_ts = size(NosepokeRew_ts,1);
A
RightNosepokeRewPsthArrayA = NaN(nRightNosepokeRew_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iRightNosepokeRew = 1:nRightNosepokeRew_ts
    thisTime = RightNosepokeRew_ts(iRightNosepokeRew)-1;
    thisIndex = round((thisTime*samplingRate));
    RightNosepokeRewPsthArrayA(iRightNosepokeRew,:) = processPhotDataRow_normDat(zDatA, thisIndex, nTsPrev, nTsPost);
end
NosepokeRewErrA1 = (nanstd(NosepokeRewPsthArrayA(1:10,1:30519)))/sqrt(size(NosepokeRewPsthArrayA(1:10,1:30519),1));
NosepokeRewPsthA1 = nanmean(NosepokeRewPsthArrayA(1:10,1:30519));
B
RightNosepokeRewPsthArrayB = NaN(nRightNosepokeRew_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iRightNosepokeRew = 1:nRightNosepokeRew_ts
    thisTime = RightNosepokeRew_ts(iRightNosepokeRew)-1;
    thisIndex = round((thisTime*samplingRate));
    RightNosepokeRewPsthArrayB(iRightNosepokeRew,:) = processPhotDataRow_normDat(zDatB, thisIndex, nTsPrev, nTsPost);
end 
NosepokeRewErrB1 = (nanstd(NosepokeRewPsthArrayB(1:10,1:30519)))/sqrt(size(NosepokeRewPsthArrayB(1:10,1:30519),1));
NosepokeRewPsthB1 = nanmean(NosepokeRewPsthArrayB(1:10,1:30519));
EarlyMeanB1 = mean(NosepokeRewPsthB1(round(10*samplingRate): round(20*samplingRate)))- mean(NosepokeRewPsthB1(round(5*samplingRate): round(10*samplingRate)));
% if plotting RightNosepokeRew
nRightNosepokeRew_ts = size(RightNosepokeRew_ts,1);
A
RightNosepokeRewPsthArrayA = NaN(nRightNosepokeRew_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iRightNosepokeRew = 1:nRightNosepokeRew_ts
    thisTime = RightNosepokeRew_ts(iRightNosepokeRew)-1;
    thisIndex = round((thisTime*samplingRate));
    RightNosepokeRewPsthArrayA(iRightNosepokeRew,:) = processPhotDataRow_normDat(zDatA, thisIndex, nTsPrev, nTsPost);
end
NosepokeRewErrA2 = (nanstd(NosepokeRewPsthArrayA(11:20,1:30519)))/sqrt(size(NosepokeRewPsthArrayA(11:20,1:30519),1));
NosepokeRewPsthA2 = nanmean(NosepokeRewPsthArrayA(11:20,1:30519));
B
RightNosepokeRewPsthArrayB = NaN(nRightNosepokeRew_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iRightNosepokeRew = 1:nRightNosepokeRew_ts
    thisTime = RightNosepokeRew_ts(iRightNosepokeRew)-1;
    thisIndex = round((thisTime*samplingRate));
    RightNosepokeRewPsthArrayB(iRightNosepokeRew,:) = processPhotDataRow_normDat(zDatB, thisIndex, nTsPrev, nTsPost);
end 
NosepokeRewErrB2 = (nanstd(NosepokeRewPsthArrayB(11:20,1:30519)))/sqrt(size(NosepokeRewPsthArrayB(11:20,1:30519),1));
NosepokeRewPsthB2 = nanmean(NosepokeRewPsthArrayB(11:20,1:30519));
EarlyMeanB2 = mean(NosepokeRewPsthB2(round(10*samplingRate): round(20*samplingRate)))- mean(NosepokeRewPsthB2(round(5*samplingRate): round(10*samplingRate)));
% if plotting RightNosepokeRew
nNosepokeRew_ts = size(NosepokeRew_ts,1);
A
RightNosepokeRewPsthArrayA = NaN(nRightNosepokeRew_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iRightNosepokeRew = 1:nRightNosepokeRew_ts
    thisTime = RightNosepokeRew_ts(iRightNosepokeRew)-1;
    thisIndex = round((thisTime*samplingRate));
    RightNosepokeRewPsthArrayA(iRightNosepokeRew,:) = processPhotDataRow_normDat(zDatA, thisIndex, nTsPrev, nTsPost);
end
NosepokeRewErrA3 = (nanstd(NosepokeRewPsthArrayA(21:30,1:30519)))/sqrt(size(NosepokeRewPsthArrayA(21:30,1:30519),1));
NosepokeRewPsthA3 = nanmean(NosepokeRewPsthArrayA(21:30,1:30519));
B
RightNosepokeRewPsthArrayB = NaN(nRightNosepokeRew_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iRightNosepokeRew = 1:nRightNosepokeRew_ts
    thisTime = RightNosepokeRew_ts(iRightNosepokeRew)-1;
    thisIndex = round((thisTime*samplingRate));
    RightNosepokeRewPsthArrayB(iRightNosepokeRew,:) = processPhotDataRow_normDat(zDatB, thisIndex, nTsPrev, nTsPost);
end 
NosepokeRewErrB3 = (nanstd(NosepokeRewPsthArrayB(21:30,1:30519)))/sqrt(size(NosepokeRewPsthArrayB(21:30,1:30519),1));
NosepokeRewPsthB3 = nanmean(NosepokeRewPsthArrayB(21:30,1:30519));
LateMeanB3 = mean(NosepokeRewPsthB3(round(10*samplingRate): round(20*samplingRate)))- mean(NosepokeRewPsthB3(round(5*samplingRate): round(10*samplingRate)));
% if plotting RightNosepokeRew
nNosepokeRew_ts = size(NosepokeRew_ts,1);
A
RightNosepokeRewPsthArrayA = NaN(nRightNosepokeRew_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iRightNosepokeRew = 1:nRightNosepokeRew_ts
    thisTime = RightNosepokeRew_ts(iRightNosepokeRew)-1;
    thisIndex = round((thisTime*samplingRate));
    RightNosepokeRewPsthArrayA(iRightNosepokeRew,:) = processPhotDataRow_normDat(zDatA, thisIndex, nTsPrev, nTsPost);
end
NosepokeRewErrA4 = (nanstd(NosepokeRewPsthArrayA(31:40,1:30519)))/sqrt(size(NosepokeRewPsthArrayA(31:40,1:30519),1));
NosepokeRewPsthA4 = nanmean(NosepokeRewPsthArrayA(31:40,1:30519));
B
RightNosepokeRewPsthArrayB = NaN(nRightNosepokeRew_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iRightNosepokeRew = 1:nRightNosepokeRew_ts
    thisTime = RightNosepokeRew_ts(iRightNosepokeRew)-1;
    thisIndex = round((thisTime*samplingRate));
    RightNosepokeRewPsthArrayB(iRightNosepokeRew,:) = processPhotDataRow_normDat(zDatB, thisIndex, nTsPrev, nTsPost);
end 
NosepokeRewErrB4 = (nanstd(NosepokeRewPsthArrayB(31:40,1:30519)))/sqrt(size(NosepokeRewPsthArrayB(31:40,1:30519),1));
NosepokeRewPsthB4 = nanmean(NosepokeRewPsthArrayB(31:40,1:30519));
LateMeanB4 = mean(NosepokeRewPsthB4(round(10*samplingRate): round(20*samplingRate)))- mean(NosepokeRewPsthB4(round(5*samplingRate): round(10*samplingRate)));
% if plotting RightNosepokeRew
nNosepokeRew_ts = size(NosepokeRew_ts,1);
A
RightNosepokeRewPsthArrayA = NaN(nRightNosepokeRew_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iRightNosepokeRew = 1:nRightNosepokeRew_ts
    thisTime = RightNosepokeRew_ts(iRightNosepokeRew)-1;
    thisIndex = round((thisTime*samplingRate));
    RightNosepokeRewPsthArrayA(iRightNosepokeRew,:) = processPhotDataRow_normDat(zDatA, thisIndex, nTsPrev, nTsPost);
end
NosepokeRewErrA5 = (nanstd(NosepokeRewPsthArrayA(41:50,1:30519)))/sqrt(size(NosepokeRewPsthArrayA(41:50,1:30519),1));
NosepokeRewPsthA5 = nanmean(NosepokeRewPsthArrayA(41:50,1:30519));
B
RightNosepokeRewPsthArrayB = NaN(nRightNosepokeRew_ts,nTsPrev+nTsPost+1); % preallocate arrays for speed
for iRightNosepokeRew = 1:nRightNosepokeRew_ts
    thisTime = RightNosepokeRew_ts(iRightNosepokeRew)-1;
    thisIndex = round((thisTime*samplingRate));
    RightNosepokeRewPsthArrayB(iRightNosepokeRew,:) = processPhotDataRow_normDat(zDatB, thisIndex, nTsPrev, nTsPost);
end 
NosepokeRewErrB5 = (nanstd(NosepokeRewPsthArrayB(41:50,1:30519)))/sqrt(size(NosepokeRewPsthArrayB(41:50,1:30519),1));
NosepokeRewPsthB5 = nanmean(NosepokeRewPsthArrayB(41:50,1:30519));
LateMeanB5 = mean(RightNosepokeRewPsthB5(round(10*samplingRate): round(20*samplingRate)))- mean(RightNosepokeRewPsthB5(round(5*samplingRate): round(10*samplingRate)));

%% make Right Nosepoke Reward plot
totalTs = nTsPrev + nTsPost;
increment = (nSecPrev + nSecPost) / totalTs;
timeAxis = (-1 * nSecPrev) : increment : nSecPost;
NosepokeRewErrPosB1 = NosepokeRewPsthB1 + NosepokeRewErrB1;
NosepokeRewErrNegB1 = NosepokeRewPsthB1 - NosepokeRewErrB1;
NosepokeRewErrPosB2 = NosepokeRewPsthB2 + NosepokeRewErrB2;
NosepokeRewErrNegB2 = NosepokeRewPsthB2 - NosepokeRewErrB2;
NosepokeRewErrPosB3 = NosepokeRewPsthB3 + NosepokeRewErrB3;
NosepokeRewErrNegB3 = NosepokeRewPsthB3 - NosepokeRewErrB3;
NosepokeRewErrPosB4 = NosepokeRewPsthB4 + NosepokeRewErrB4;
NosepokeRewErrNegB4 = NosepokeRewPsthB4 - NosepokeRewErrB4;
NosepokeRewErrPosB5 = NosepokeRewPsthB5 + NosepokeRewErrB5;
NosepokeRewErrNegB5 = NosepokeRewPsthB5 - NosepokeRewErrB5;

figure('Name', 'Right Nosepoke Reward', 'Color',[1 1 1])
hold on
set(gca, 'Fontsize', 18);
xlim ([-5,10]);
ylim ([-2, 4]);
h1 = plot(timeAxis, NosepokeRewPsthB1, 'r-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[NosepokeRewErrPosB1, fliplr(NosepokeRewErrNegB1)], [1,0,0],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis, NosepokeRewPsthB2, 'color', [.96, .65, .26], 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[NosepokeRewErrPosB2, fliplr(NosepokeRewErrNegB2)], [.96, .65, .26],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h3 = plot(timeAxis, NosepokeRewPsthB3, 'color', [.23, .66,.16],  'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[NosepokeRewErrPosB3, fliplr(NosepokeRewErrNegB3)], [.23, .66,.16],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h4 = plot(timeAxis, NosepokeRewPsthB4, 'b-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[NosepokeRewErrPosB4, fliplr(NosepokeRewErrNegB4)], [0,0,1],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h5 = plot(timeAxis, NosepokeRewPsthB5, 'color', [.46, .09, .69], 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[NosepokeRewErrPosB5, fliplr(NosepokeRewErrNegB5)], [.46, .09, .69],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontSize', 18);
ylabel('Z-score', 'FontSize', 18);
legend ([h1, h2], '1-10', '11-21', 'orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;
