%% Extract photometry data from TDT Data Tank
%UPDATED Oct 2017 for Lerner Lab

%path_to_data = 'R:\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\Lerner_Lab_tnl2633\Talia\FiberPhotometry'; % point to folder where data is
path_to_data = 'C:\Users\tnl2633\Documents';

% set up Tank name, variables to extract
tankdir = path_to_data;
tankname = 'DEMOTANK2';
blockname = 'test10262017';
storenames = {'Dv1A' 'Dv2A' 'Dv3B' 'Dv4B' 'Ep1/' 'Ep2/' 'Ep3/'}; % name of stores to extract from TDT (usu. 4-letter code) LMag is the demodulated data, may also have other timestamps etc
% extract
for k = 1:numel(storenames)
  storename = storenames{k};
  S{k} = tdt2mat(tankdir, tankname, blockname, storename);
end

%% Massage data and get time stamps

%add more if you extracted more stores above
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
