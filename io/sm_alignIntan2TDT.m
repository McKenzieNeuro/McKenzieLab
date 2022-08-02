function output = sm_alignIntan2TDT(filepath)
fs = 30000; %read from Intan xml


fils = dir(filepath);
kp_TDT = cellfun(@any,regexp({fils.name},'-'));
kp_Intan = cellfun(@any,regexp({fils.name},'_'));

IntanDir = fils(kp_Intan);
TDTDir = fils(kp_TDT);
%
data = TDTbin2mat(TDTDir);
ts = (1:length(data.streams.Dv4D.data)) / data.streams.Dv4D.fs;
TDT_pulsein = data.epocs.PC0_.onset;
digitalIn = [filepath filesep IntanDir.name filesep 'digitalin.dat'];
digIn = LoadBinary(digitalIn,'nchannels',1);
pulses_intan = find(diff(digIn)>0);

kp = diff(pulses_intan)/fs > .99 & diff(pulses_intan)/fs < 1.01;

pulses_intan = pulses_intan(kp);

if length(TDT_pulsein)<length(pulses_intan)
    %grab just the last pulses, likely pulses were sent during preview
    pulses_intan = pulses_intan(end-length(TDT_pulsein)+1:end);
    
end

pulses_intan = pulses_intan/fs;

output = [pulses_intan TDT_pulsein];

