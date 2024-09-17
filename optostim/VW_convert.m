function V_out = VW_convert(maxstimV,dur,shank,LED,varargin)

%define half sine stimulation of a particular duration
warning off

VWfname = [];
outfil = 'halfsine.txt';
dt = 1/1000;
minV = 2;
% Parse options
for i = 1:2:length(varargin),
    if ~isa(varargin{i},'char'),
        error(['Parameter ' num2str(i+3) ' is not a property ']);
    end
    switch(lower(varargin{i})),
        case 'vwfname',
            VWfname = varargin{i+1};
            if ~isa(VWfname,'char')
                error('Incorrect value for property ''VWfname''');
            end
        case 'dt'
            dt = varargin{i+1};
            if ~isa(dt,'double')
                error('Incorrect value for property ''dt''');
            end
        case 'outfil',
            outfil = varargin{i+1};
            if ~isa(outfil,'char')
                error('Incorrect value for property ''VWfname''');
            end
            
        case 'minV',
            minV = varargin{i+1};
            if ~isa(minV,'double')
                error('Incorrect value for property ''VWfname''');
            end
            
            
            
            
    end
end



if isempty(VWfname)
    [VWfname,dirName] = uigetfile('/home/sam/Dropbox/buzsaki/PROTOCOLS/OptData/calibration_files/*.xls','Choose VW conversion file');
else
    dirName = pwd;
end



%% Conversion to the current prep
WV = xlsread([dirName VWfname],['S' num2str(shank) 'L' num2str(LED)],'A1:B6');
V = WV(:, 1);
W = [WV(:, 2)];


V_int = [minV:.001:max(V)];
W_int = interp1([minV;V],[0;W],V_int,'cubic');

%%
%get max of half sine

[~,ind] = bestmatch(maxstimV,V_int);
maxstimW = W_int(ind);


%define output light waveform
k =  maxstimW/2 * -cos([dt:dt:dur-dt]*2*pi/dur)+maxstimW/2;
ind = nan(length(ind)-1,1);
for i = 1:length(k)
    [~,ind(i)] = bestmatch(k(i),W_int);
end

V_out = V_int(ind);

V_out(V_out<2.5) = [];

    entrance = exp(  log(.001):.2:log(.5))+2;




V_out = [entrance V_out fliplr(entrance)];
%%
outfil = [dirName '/' outfil];
dlmwrite(outfil,V_out(:))

end

