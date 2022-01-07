%% Combine photometry traces from different mice

%Copy the relevant .mat files together into a folder

% basename = 'LNH';
folderName=uigetdir('C:\');
d = dir(folderName);

f_names={};
for i=1:length(d)
    fn=d(i).name;
    f_names=[f_names; fn];
end
mat_idx= find(contains(f_names(:,1), '.mat'));
f_names=f_names(mat_idx);
Number_mat=length(f_names);

% look for .mat files            % find the number of .mat files to run in the for-loop   

load([f_names{1}], 'nTsPrev', 'nTsPost');     % load in PSTH window variables from the first .mat file to use for array preallocations (assuming these are same for all files)
%%
% Decide which vectors you want to extract and preallocate arrays for speed
% Preallocate Right Reward
GroupNosepokeRewPsthArrayA = NaN (Number_mat, nTsPrev+nTsPost+1);
GroupNosepokeRewPsthArrayB = NaN (Number_mat, nTsPrev+nTsPost+1);

%% 
%Preallocate Right No Reward
GroupNosepokeNoRewPsthArrayA = NaN (Number_mat, nTsPrev+nTsPost+1);
GroupNosepokeNoRewPsthArrayB = NaN (Number_mat, nTsPrev+nTsPost+1);

%%
%Preallocate Left Reward
GroupControlPsthArrayA = NaN (Number_mat, nTsPrev+nTsPost+1);
GroupControlPsthArrayB = NaN (Number_mat, nTsPrev+nTsPost+1);

%%
%Preallocate Port Entry Reward
GroupPortRewardPsthArrayA = NaN (Number_mat, nTsPrev+nTsPost+1);
GroupPortRewardPsthArrayB = NaN (Number_mat, nTsPrev+nTsPost+1);

%%
%Preallocate Port Entry No Reward
GroupPortNoRewardPsthArrayA = NaN (Number_mat, nTsPrev+nTsPost+1);
GroupPortNoRewardPsthArrayB = NaN (Number_mat, nTsPrev+nTsPost+1);

%%
GroupShockPsthArrayA = NaN (Number_mat, nTsPrev+nTsPost+1);
GroupShockPsthArrayB = NaN (Number_mat, nTsPrev+nTsPost+1);

%%
% Get data from each .mat file to fill arrays

for i = 1: length(f_names)
load(f_names{i});
 GroupNosepokeRewPsthArrayA (i,:) = NosepokeRewPsthA; %put the PSTHs from the file into row i of the new group arrays
 GroupNosepokeRewPsthArrayB (i,:) = NosepokeRewPsthB;
end

%%
% Get Right No Reward data from each .mat file
for i = 1: length(f_names)
 load(f_names{i});%load a .mat file
 GroupNosepokeNoRewPsthArrayA (i,:) = NosepokeNoRewPsthA; %put the PSTHs from the file into row i of the new group arrays
 GroupNosepokeNoRewPsthArrayB (i,:) = NosepokeNoRewPsthB;
end

%%
% Get Left Reward data from each .mat file
for i = 1: length(f_names)
 load(f_names{i});
 GroupControlPsthArrayA (i,:) = ControlPsthA; %put the PSTHs from the file into row i of the new group arrays
 GroupControlPsthArrayB (i,:) = ControlPsthB;
end

%%
% Get Port Reward data from each .mat file
for i = 1: length(f_names)
 %load([basename '_' num2str(i) '.mat'])                         %load a .mat file
 load(f_names{i});
 GroupPortRewardPsthArrayA (i,:) = PortRewardPsthA; %put the PSTHs from the file into row i of the new group arrays
 GroupPortRewardPsthArrayB (i,:) = PortRewardPsthB;
end

%%
% Get Port No Reward data from each .mat file
for i = 1: length(f_names)
 %load([basename '_' num2str(i) '.mat']) 
 load(f_names{i}); %load a .mat file
 GroupPortNoRewardPsthArrayA (i,:) = PortNoRewardPsthA; %put the PSTHs from the file into row i of the new group arrays
 GroupPortNoRewardPsthArrayB (i,:) = PortNoRewardPsthB;
end

%%
% Get Shock data from each .mat file
for i = 1: length(f_names)
 %load([basename '_' num2str(i) '.mat'])                         %load a .mat file
 load(f_names{i});
 GroupShockPsthArrayA (i,:) = ShockPsthA; %put the PSTHs from the file into row i of the new group arrays
 GroupShockPsthArrayB (i,:) = ShockPsthB;
end


%%
% Combine Right Reward data into one vector
GroupNosepokeRewPsthA = nanmean (GroupNosepokeRewPsthArrayA);
GroupNosepokeRewErrA = (nanstd(GroupNosepokeRewPsthArrayA))/sqrt(size(GroupNosepokeRewPsthArrayA,1));
GroupNosepokeRewPsthB = nanmean (GroupNosepokeRewPsthArrayB);
GroupNosepokeRewErrB = (nanstd(GroupNosepokeRewPsthArrayB))/sqrt(size(GroupNosepokeRewPsthArrayB,1));

%... etc for each thing you're combining

%%
% Combine Right No Reward data into one vector
GroupNosepokeNoRewPsthA = nanmean (GroupNosepokeNoRewPsthArrayA);
GroupNosepokeNoRewErrA = (nanstd(GroupNosepokeNoRewPsthArrayA))/sqrt(size(GroupNosepokeNoRewPsthArrayA,1));
GroupNosepokeNoRewPsthB = nanmean (GroupNosepokeNoRewPsthArrayB);
GroupNosepokeNoRewErrB = (nanstd(GroupNosepokeNoRewPsthArrayB))/sqrt(size(GroupNosepokeNoRewPsthArrayB,1));
%% 
% Combine Left Reward data into one vector
GroupControlPsthA = nanmean (GroupControlPsthArrayA);
GroupControlErrA = (nanstd(GroupControlPsthArrayA))/sqrt(size(GroupControlPsthArrayA,1));
GroupControlPsthB = nanmean (GroupControlPsthArrayB);
GroupControlErrB = (nanstd(GroupControlPsthArrayB))/sqrt(size(GroupControlPsthArrayB,1));

%%
% Combine Port Reward data into one vector
GroupPortRewardPsthA = nanmean (GroupPortRewardPsthArrayA);
GroupPortRewardErrA = (nanstd(GroupPortRewardPsthArrayA))/sqrt(size(GroupPortRewardPsthArrayA,1));
GroupPortRewardPsthB = nanmean (GroupPortRewardPsthArrayB);
GroupPortRewardErrB = (nanstd(GroupPortRewardPsthArrayB))/sqrt(size(GroupPortRewardPsthArrayB,1));
%%
% Combine Port No Reward data into one vector
GroupPortNoRewardPsthA = nanmean (GroupPortNoRewardPsthArrayA);
GroupPortNoRewardErrA = (nanstd(GroupPortNoRewardPsthArrayA))/sqrt(size(GroupPortNoRewardPsthArrayA,1));
GroupPortNoRewardPsthB = nanmean (GroupPortNoRewardPsthArrayB);
GroupPortNoRewardErrB = (nanstd(GroupPortNoRewardPsthArrayB))/sqrt(size(GroupPortNoRewardPsthArrayB,1));

%%
% Combine Shock data into one vector
GroupShockPsthA = nanmean (GroupShockPsthArrayA);
GroupShockErrA = (nanstd(GroupShockPsthArrayA))/sqrt(size(GroupShockPsthArrayA,1));
GroupShockPsthB = nanmean (GroupShockPsthArrayB);
GroupShockErrB = (nanstd(GroupShockPsthArrayB))/sqrt(size(GroupShockPsthArrayB,1));

%%
%Plot
%make Rewarded Right Nosepoke plot for combined data
totalTs = nTsPrev + nTsPost;
increment = (nSecPrev + nSecPost) / totalTs;
timeAxis = (-1 * nSecPrev) : increment : nSecPost;
GroupNosepokeRewErrPosA = GroupNosepokeRewPsthA + GroupNosepokeRewErrA;
GroupNosepokeRewErrNegA = GroupNosepokeRewPsthA - GroupNosepokeRewErrA;
GroupNosepokeRewErrPosB = GroupNosepokeRewPsthB + GroupNosepokeRewErrB;
GroupNosepokeRewErrNegB = GroupNosepokeRewPsthB - GroupNosepokeRewErrB;
figure('Name', 'Rewarded Nosepokes','Color',[1 1 1])
hold on
set(gca, 'FontName','Arial','FontSize',26,'FontSmoothing','off',...
    'FontWeight','bold','LineWidth',4);
xlim ([-5,10]);
ylim ([-0.5,3.5]);
h1 = plot(timeAxis,GroupNosepokeRewPsthA, 'r-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupNosepokeRewErrPosA, fliplr(GroupNosepokeRewErrNegA)], [1,0,0],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis,GroupNosepokeRewPsthB, 'b-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupNosepokeRewErrPosB, fliplr(GroupNosepokeRewErrNegB)], [0,0,1],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontName','Arial','FontSize',26,'FontSmoothing','off',...
    'FontWeight','bold','LineWidth',4);
ylabel('z-score', 'FontName','Arial','FontSize',26,'FontSmoothing','off',...
    'FontWeight','bold','LineWidth',4);
legend ([h1, h2], 'DMS','DLS','orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;
print ('RewardedNosepokesHR', '-dsvg', '-r300');
% 
%%
%make No Reward Right Nosepoke plot for combined data
totalTs = nTsPrev + nTsPost;
increment = (nSecPrev + nSecPost) / totalTs;
timeAxis = (-1 * nSecPrev) : increment : nSecPost;
GroupNosepokeNoRewErrPosA = GroupNosepokeNoRewPsthA + GroupNosepokeNoRewErrA;
GroupNosepokeNoRewErrNegA = GroupNosepokeNoRewPsthA - GroupNosepokeNoRewErrA;
GroupNosepokeNoRewErrPosB = GroupNosepokeNoRewPsthB + GroupNosepokeNoRewErrB;
GroupNosepokeNoRewErrNegB = GroupNosepokeNoRewPsthB - GroupNosepokeNoRewErrB;
figure('Name', 'Unrewarded Nosepokes','Color',[1 1 1]);
hold on
set(gca, 'Fontsize', 14);
xlim ([-5,10]);
h1 = plot(timeAxis,GroupNosepokeNoRewPsthA, 'r-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupNosepokeNoRewErrPosA, fliplr(GroupNosepokeNoRewErrNegA)], [1,0,0],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis,GroupNosepokeNoRewPsthB, 'b-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupNosepokeNoRewErrPosB, fliplr(GroupNosepokeNoRewErrNegB)], [0,0,1],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontSize', 12);
ylabel('z-score', 'FontSize', 12);
legend ([h1, h2], 'DMS','DLS','orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;
print ('UnrewardedNosepokesHR', '-dsvg', '-r300');

% 
 %%
%make Control plot for combined data
totalTs = nTsPrev + nTsPost;
increment = (nSecPrev + nSecPost) / totalTs;
timeAxis = (-1 * nSecPrev) : increment : nSecPost;
GroupControlErrPosA = GroupControlPsthA + GroupControlErrA;
GroupControlErrNegA = GroupControlPsthA - GroupControlErrA;
GroupControlErrPosB = GroupControlPsthB + GroupControlErrB;
GroupControlErrNegB = GroupControlPsthB - GroupControlErrB;
figure('Name', 'Control Nosepokes','Color',[1 1 1])
hold on
set(gca, 'Fontsize', 14);
xlim ([-5,10]);
h1 = plot(timeAxis,GroupControlPsthA, 'color', [.95 .45 .17], 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupControlErrPosA, fliplr(GroupControlErrNegA)], [0,0,1],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis,GroupControlPsthB, 'color', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupControlErrPosB, fliplr(GroupControlErrNegB)], [1,0,0],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontSize', 12);
ylabel('z-score', 'FontSize', 12);
legend ([h1, h2], 'DMS','DLS','orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;
% 
%%
%make Rewarded Port Entry plot for combined data
totalTs = nTsPrev + nTsPost;
increment = (nSecPrev + nSecPost) / totalTs;
timeAxis = (-1 * nSecPrev) : increment : nSecPost;
GroupPortRewardErrPosA = GroupPortRewardPsthA + GroupPortRewardErrA;
GroupPortRewardErrNegA = GroupPortRewardPsthA - GroupPortRewardErrA;
GroupPortRewardErrPosB = GroupPortRewardPsthB + GroupPortRewardErrB;
GroupPortRewardErrNegB = GroupPortRewardPsthB - GroupPortRewardErrB;
figure('Name', 'Rewarded Port Entry','Color',[1 1 1])
hold on
set(gca, 'FontName','Arial','FontSize',26,'FontSmoothing','off',...
    'FontWeight','bold','LineWidth',4);
xlim ([-5,10]);
ylim ([-.5,2.5])
h1 = plot(timeAxis,GroupPortRewardPsthA, 'color', [.55 .04 .07], 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupPortRewardErrPosA, fliplr(GroupPortRewardErrNegA)], [.55 .04 .07],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis,GroupPortRewardPsthB, 'color', [0 .18 .28], 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupPortRewardErrPosB, fliplr(GroupPortRewardErrNegB)], [0 .18 .28],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontName','Arial','FontSize',26,'FontSmoothing','off',...
    'FontWeight','bold','LineWidth',4);
ylabel('z-score', 'FontName','Arial','FontSize',26,'FontSmoothing','off',...
    'FontWeight','bold','LineWidth',4);
legend ([h1, h2], 'DMS','DLS','orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;
print ('RewardedPortHR', '-dsvg', '-r300');

%%
%make Unrewarded Port Entry plot for combined data
totalTs = nTsPrev + nTsPost;
increment = (nSecPrev + nSecPost) / totalTs;
timeAxis = (-1 * nSecPrev) : increment : nSecPost;
GroupPortNoRewardErrPosA = GroupPortNoRewardPsthA + GroupPortNoRewardErrA;
GroupPortNoRewardErrNegA = GroupPortNoRewardPsthA - GroupPortNoRewardErrA;
GroupPortNoRewardErrPosB = GroupPortNoRewardPsthB + GroupPortNoRewardErrB;
GroupPortNoRewardErrNegB = GroupPortNoRewardPsthB - GroupPortNoRewardErrB;
figure('Name', 'Unewarded Port Entry','Color',[1 1 1])
hold on
set(gca, 'FontName','Arial','FontSize',26,'FontSmoothing','off',...
    'FontWeight','bold','LineWidth',4);
xlim ([-5,10]);
ylim ([-.75,.75]);
h1 = plot(timeAxis,GroupPortNoRewardPsthA, 'color', [.95 .45 .17], 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupPortNoRewardErrPosA, fliplr(GroupPortNoRewardErrNegA)], [.95 .45 .17],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis,GroupPortNoRewardPsthB, 'color', [.24 .73 .92], 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupPortNoRewardErrPosB, fliplr(GroupPortNoRewardErrNegB)], [.24 .73 .92],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontName','Arial','FontSize',26,'FontSmoothing','off',...
    'FontWeight','bold','LineWidth',4);
ylabel('z-score', 'FontName','Arial','FontSize',26,'FontSmoothing','off',...
    'FontWeight','bold','LineWidth',4);
legend ([h1, h2], 'DMS','DLS','orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;
print ('UnrewardedPortHR', '-dsvg', '-r300');

%%
%make Shock plot for combined data
totalTs = nTsPrev + nTsPost;
increment = (nSecPrev + nSecPost) / totalTs;
timeAxis = (-1 * nSecPrev) : increment : nSecPost;
GroupShockErrPosA = GroupShockPsthA + GroupShockErrA;
GroupShockErrNegA = GroupShockPsthA - GroupShockErrA;
GroupShockErrPosB = GroupShockPsthB + GroupShockErrB;
GroupShockErrNegB = GroupShockPsthB - GroupShockErrB;
figure('Name', 'Shocks','Color',[1 1 1])
hold on
set(gca, 'Fontsize', 14);
xlim ([-5,10]);
h1 = plot(timeAxis,GroupShockPsthA, 'b-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupShockErrPosA, fliplr(GroupShockErrNegA)], [0,0,1],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis,GroupShockPsthB, 'r-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupShockErrPosB, fliplr(GroupShockErrNegB)], [1,0,0],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontSize', 12);
ylabel('z-score', 'FontSize', 12);
legend ([h1, h2], 'DMS','DLS','orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;
% %%
totalTs = nTsPrev + nTsPost;
increment = (nSecPrev + nSecPost) / totalTs;
timeAxis = (-1 * nSecPrev) : increment : nSecPost;
GroupReward1ErrPosB = GroupReward1PsthB + GroupReward1ErrB;
GroupReward1ErrNegA = GroupReward1PsthB - GroupReward1ErrB;
GroupReward2ErrPosA = GroupReward2PsthB + GroupReward2ErrB;
GroupReward2ErrNegA = GroupReward2PsthB - GroupReward2ErrB;
GroupReward3ErrPosA = GroupReward3PsthB + GroupReward3ErrB;
GroupReward3ErrNegA = GroupReward3PsthB - GroupReward3ErrB;
GroupReward4ErrPosA = GroupReward4PsthB + GroupReward4ErrB;
GroupReward4ErrNegA = GroupReward4PsthB - GroupReward4ErrB;
GroupReward5ErrPosA = GroupReward5PsthB + GroupReward5ErrB;
GroupReward5ErrNegA = GroupReward5PsthB - GroupReward5ErrB;

figure('Name', 'Right Nosepoke Reward', 'Color',[1 1 1])
hold on
set(gca, 'Fontsize', 18);
xlim ([-5,10]);
ylim([0,.5]);
h1 = plot(timeAxis, GroupReward1PsthB, 'r-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupReward1ErrPosB, fliplr(GroupReward1ErrNegA)], [1,0,0],'FaceAlpha', 0.2, 'EdgeColor', 'none');
% h2 = plot(timeAxis, GroupReward2PsthB, 'color', [.96, .65, .26], 'linewidth', 2);
% fill([timeAxis, fliplr(timeAxis)],[GroupReward2ErrPosA, fliplr(GroupReward2ErrNegA)], [.96, .65, .26],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis, GroupReward3PsthB, 'color', [.23, .66,.16],  'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupReward3ErrPosA, fliplr(GroupReward3ErrNegA)], [.23, .66,.16],'FaceAlpha', 0.2, 'EdgeColor', 'none');
% h4 = plot(timeAxis, GroupReward4PsthB, 'b-', 'linewidth', 2);
% fill([timeAxis, fliplr(timeAxis)],[GroupReward4ErrPosA, fliplr(GroupReward4ErrNegA)], [0,0,1],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h3 = plot(timeAxis, GroupReward5PsthB, 'color', [.46, .09, .69], 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupReward5ErrPosA, fliplr(GroupReward5ErrNegA)], [.46, .09, .69],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontSize', 18);
ylabel('Z-score', 'FontSize', 18);
legend ([h1, h2], '1-10', '11-21', 'orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;

%%
totalTs = nTsPrev + nTsPost;
increment = (nSecPrev + nSecPost) / totalTs;
timeAxis = (-1 * nSecPrev) : increment : nSecPost;
GroupReward1ErrPosA = GroupReward1PsthA + GroupReward1ErrA;
GroupReward1ErrNegA = GroupReward1PsthA - GroupReward1ErrA;
GroupReward2ErrPosA = GroupReward2PsthA + GroupReward2ErrA;
GroupReward2ErrNegA = GroupReward2PsthA - GroupReward2ErrA;
GroupReward3ErrPosA = GroupReward3PsthA + GroupReward3ErrA;
GroupReward3ErrNegA = GroupReward3PsthA - GroupReward3ErrA;
GroupReward4ErrPosA = GroupReward4PsthA + GroupReward4ErrA;
GroupReward4ErrNegA = GroupReward4PsthA - GroupReward4ErrA;
GroupReward5ErrPosA = GroupReward5PsthA + GroupReward5ErrA;
GroupReward5ErrNegA = GroupReward5PsthA - GroupReward5ErrA;

figure('Name', 'Right Nosepoke Reward', 'Color',[1 1 1])
hold on
set(gca, 'Fontsize', 18);
xlim ([-5,10]);
h1 = plot(timeAxis, GroupReward1PsthA, 'r-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupReward1ErrPosA, fliplr(GroupReward1ErrNegA)], [1,0,0],'FaceAlpha', 0.2, 'EdgeColor', 'none');
% h2 = plot(timeAxis, GroupReward2PsthA, 'color', [.96, .65, .26], 'linewidth', 2);
% fill([timeAxis, fliplr(timeAxis)],[GroupReward2ErrPosA, fliplr(GroupReward2ErrNegA)], [.96, .65, .26],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis, GroupReward3PsthA, 'color', [.23, .66,.16],  'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupReward3ErrPosA, fliplr(GroupReward3ErrNegA)], [.23, .66,.16],'FaceAlpha', 0.2, 'EdgeColor', 'none');
% h4 = plot(timeAxis, GroupReward4PsthA, 'b-', 'linewidth', 2);
% fill([timeAxis, fliplr(timeAxis)],[GroupReward4ErrPosA, fliplr(GroupReward4ErrNegA)], [0,0,1],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h5 = plot(timeAxis, GroupReward5PsthA, 'color', [.46, .09, .69], 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupReward5ErrPosA, fliplr(GroupReward5ErrNegA)], [.46, .09, .69],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontSize', 18);
ylabel('Z-score', 'FontSize', 18);
legend ([h1, h2], '1-10', '11-21', 'orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;

%%
figure('Name', 'Right Nosepoke Reward', 'Color',[1 1 1])
hold on
set(gca, 'FontName','Arial','FontSize',26,'FontSmoothing','off',...
    'FontWeight','bold','LineWidth',4);
xlim ([-5,10]);
ylim ([-0.6,4]);
h1 = plot(timeAxis,GroupNosepokeRewPsthA, 'color', [.55 .04 .07],'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupNosepokeRewErrPosA, fliplr(GroupNosepokeRewErrNegA)], [.55 .04 .07],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis,GroupNosepokeNoRewPsthA, 'color', [.95 .45 .17], 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupNosepokeNoRewErrPosA, fliplr(GroupNosepokeNoRewErrNegA)], [.95 .45 .17],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontName','Arial','FontSize',26,'FontSmoothing','off',...
    'FontWeight','bold','LineWidth',4);
ylabel('z-score', 'FontName','Arial','FontSize',26,'FontSmoothing','off',...
    'FontWeight','bold','LineWidth',4);
legend ([h1, h2], 'Rewarded','Unrewarded','orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;
print ('DMSNosepokesHR', '-dsvg', '-r300');

%%
figure('Name', 'Right Nosepoke Reward', 'Color',[1 1 1])
hold on
set(gca, 'FontName','Arial','FontSize',26,'FontSmoothing','off',...
    'FontWeight','bold','LineWidth',4);
xlim ([-5,10]);
ylim ([-.8, 2]);
h1 = plot(timeAxis,GroupNosepokeRewPsthA, 'color', [0 .18 .28],'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupNosepokeRewErrPosA, fliplr(GroupNosepokeRewErrNegA)], [0 .18 .28],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis,GroupNosepokeNoRewPsthA, 'color', [.24 .73 .92], 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupNosepokeNoRewErrPosA, fliplr(GroupNosepokeNoRewErrNegA)], [.24 .73 .92],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h3 = plot(timeAxis,GroupPortRewardPsthB, 'color', [.55 .04 .07], 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupPortRewardErrPosB, fliplr(GroupPortRewardErrNegB)], [.55 .04 .07],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h4 = plot(timeAxis,GroupPortNoRewardPsthB, 'color', [.95 .45 .17], 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupPortNoRewardErrPosA, fliplr(GroupPortNoRewardErrNegB)], [.95 .45 .17],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontName','Arial','FontSize',26,'FontSmoothing','off',...
    'FontWeight','bold','LineWidth',4);
ylabel('z-score', 'FontName','Arial','FontSize',26,'FontSmoothing','off',...
    'FontWeight','bold','LineWidth',4);
legend ([h1, h2], 'Rewarded','Unrewarded','orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;
print ('-painters','NosepokeHR', '-dsvg', '-r300');
print ('-painters','NosepokeHR', '-dtiff', '-r300');

%%
figure('Name', 'Rewarded Port Entry','Color',[1 1 1])
hold on
set(gca, 'Fontsize', 26);
xlim ([-5,10]);
ylim ([-.6, 2.2]);
h1 = plot(timeAxis,GroupPortRewardPsthA, 'color', [.55 .04 .07], 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupPortRewardErrPosA, fliplr(GroupPortRewardErrNegA)], [.55 .04 .07],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis,GroupPortNoRewardPsthA, 'color', [.95 .45 .17], 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupPortNoRewardErrPosA, fliplr(GroupPortNoRewardErrNegA)], [.95 .45 .17],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontSize', 26);
ylabel('z-score', 'FontSize', 26);
legend ([h1, h2], 'Rewarded','Unrewarded','orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;
print ('DMSPortHR', '-dsvg', '-r300');
% 

%%
figure('Name', 'Rewarded Port Entry','Color',[1 1 1])
hold on
set(gca, 'Fontsize', 26);
xlim ([-5,10]);
ylim ([-.6, 1.5]);
h1 = plot(timeAxis,GroupPortRewardPsthB, 'color', [0 .18 .28], 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupPortRewardErrPosB, fliplr(GroupPortRewardErrNegB)], [0 .18 .28],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis,GroupPortNoRewardPsthB, 'color', [.24 .73 .92], 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupPortNoRewardErrPosB, fliplr(GroupPortNoRewardErrNegB)], [.24 .73 .92],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontSize', 26);
ylabel('z-score', 'FontSize', 26);
legend ([h1, h2], 'Rewarded','Unrewarded', 'orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;
print ('DLSPortHR', '-dsvg', '-r300');

%% make Right Nosepoke Reward plot
totalTs = nTsPrev + nTsPost;
increment = (nSecPrev + nSecPost) / totalTs;
timeAxis = (-1 * nSecPrev) : increment : nSecPost;
GroupNosepokeRewErrPosA = GroupNosepokeRewPsthA + GroupNosepokeRewErrA;
GroupNosepokeRewErrNegA = GroupNosepokeRewPsthA - GroupNosepokeRewErrA;
GroupNosepokeRewErrPosB = GroupNosepokeRewPsthB + GroupNosepokeRewErrB;
GroupNosepokeRewErrNegB = GroupNosepokeRewPsthB - GroupNosepokeRewErrB;
GroupNosepokeRewPLOTPsthA = GroupNosepokeRewPsthA - GroupShockPsthA;
GroupNosepokeRewPLOTPsthB = GroupNosepokeRewPsthB - GroupShockPsthB;
GroupShockErrPosA = GroupShockPsthA + GroupShockErrA;
GroupShockErrNegA = GroupShockPsthA - GroupShockErrA;
GroupShockErrPosB = GroupShockPsthB + GroupShockErrB;
GroupShockErrNegB = GroupShockPsthB - GroupShockErrB;
figure('Name', 'Right Nosepoke Reward', 'Color',[1 1 1])
hold on
set(gca, 'Fontsize', 14);
xlim ([-5,10]);
h1 = plot(timeAxis,GroupNosepokeRewPLOTPsthA, 'b-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[(GroupNosepokeRewErrPosA - GroupShockErrPosA), fliplr(GroupNosepokeRewErrNegA - GroupShockErrNegA)], [0,0,1],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis,GroupNosepokeRewPLOTPsthB, 'r-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[(GroupNosepokeRewErrPosB - GroupShockErrPosB), fliplr(GroupNosepokeRewErrNegB- GroupShockErrNegB)], [1,0,0],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontSize', 12);
ylabel('Z-score', 'FontSize', 12);
legend ([h1, h2], 'Sensor A','Sensor B', 'orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;

%%
figure('Name', 'Shocks','Color',[1 1 1])
hold on
set(gca, 'Fontsize', 14);
xlim ([-5,10]);
h1 = plot(timeAxis,GroupShockPsthA, 'b-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupShockErrPosA, fliplr(GroupShockErrNegA)], [0,0,1],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis,GroupNosepokeRewPLOTPsthA, 'color', [0 .18 .28], 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[(GroupNosepokeRewErrPosA - GroupShockErrPosA), fliplr(GroupNosepokeRewErrNegA - GroupShockErrNegA)], [0 .18 .28],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontSize', 12);
ylabel('Z-score', 'FontSize', 12);
legend ([h1, h2], 'Shocks','Nosepoke Rewards','orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;

%%
figure('Name', 'Shocks','Color',[1 1 1])
hold on
set(gca, 'Fontsize', 14);
xlim ([-5,10]);
h1 = plot(timeAxis,GroupShockPsthB, 'b-', 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[GroupShockErrPosB, fliplr(GroupShockErrNegB)], [0,0,1],'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(timeAxis,GroupNosepokeRewPLOTPsthB, 'color',  [0 .18 .28], 'linewidth', 2);
fill([timeAxis, fliplr(timeAxis)],[(GroupNosepokeRewErrPosB - GroupShockErrPosB), fliplr(GroupNosepokeRewErrNegB - ShockErrNegB)], [0 .18 .28],'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)', 'FontSize', 12);
ylabel('Z-score', 'FontSize', 12);
legend ([h1, h2], 'Shocks','Nosepoke Rewards','orientation', 'vertical', 'Location', 'NorthEast');
legend BOXOFF;