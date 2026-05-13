%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Seizure Detection Model - Troubleshooting
            % Last Modified 26 Jan 2026 by KLuchini

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sam's Code - 1st Canoncial file
% 
cd('R:\DGregg\NeuralData\PCP\Recordings\BaseLine\8-8-2024(17.31)\RHS_240808_173241') % change current folder

v = load('auto_sz_classifier.mat');  % load model
v1 = load('PCP 4.0.sessiondata.mat'); % load ground truth

ts  = 1:size(v.estimateLabel,1);    % time vector

sz_ep = [cell2mat({v1.sessiondata.events.start}') cell2mat({v1.sessiondata.events.stop}')]; % load Ground Truth start and end pts for confirmed seizure events

kp = InIntervals(ts,sz_ep);     % kp == GT sz indices
    % tests which values fall in a list of intervals (start/stop pairs) == logical response
% 
% auto_sz =  v.estimateLabel(:,1)==2; % events detected and classified by model as sz
%     % class 2 = seizure
% 
% ons = find(diff([0;auto_sz])>0);    % find detected sz "ons"
% offs = find(diff([0;auto_sz])<0);   % find detected sz "offs"
% 
% 
% % if recording stopped mid-seizure, offs = end of recording
% if length(ons) > length(offs)
%     maxT = sm_getFileDur(v1.sessiondata.lfp_file); % from GT data (v1)
%     offs = [offs;maxT];
% end
% 
% % find detected sz events >10s
% sz_ep_auto = [ons offs];            
% kp_auto = diff(sz_ep_auto,[],2)>10;
% 
% sz_ep_auto = sz_ep_auto(kp_auto,:);
% kp_auto = InIntervals(ts,sz_ep_auto);
% 
% %detected as seizure
% seiz_auto_short = v.estimateLabel(:,1)==2 & ~kp_auto; % logical operation == detected sz AND NOT detected sz >10s
% 
% scores = v.estimateLabel(:,3); % col 3 =  prob from cat 2 = seizure
% 
% scores = scores(:);   % ensure column vectors
% scoresDur = scores;   % ensure column vectors
%     % ASK ALEX WHAT PURPOSE THE ABOVE FUNCTIONS SERVE
%       % why have scoresDur = scores?
% 
% 
% scoresDur(seiz_auto_short) = 0; % hard code scores of short auto seizures to 0
% 
% labels = double(kp); % type cast binary into double
% 
% 
% % Basic validation
% assert(numel(scores) == numel(labels), ...
%     'Scores and labels must be the same length.');

%% Compute ROC curve

% X = false positive rate
% Y = true positive rate
% T = decision thresholds
% AUC = area under the ROC curve

[Xdur, Ydur, Tdur, AUCdur] = perfcurve(labels, scoresDur, 1);
[X, Y, T, AUC] = perfcurve(labels, scores, 1);

% Plot ROC
figure;
plot(X, Y, 'LineWidth', 2);
hold on;
plot(Xdur, Ydur, 'LineWidth', 2);
plot([0 1], [0 1], 'k--');   % chance line
hold off;

xlabel('False Positive Rate');
ylabel('True Positive Rate');
title(sprintf('ROC Curve (AUC = %.3f)', AUC));
grid on;

%% Kestrel's Attempt at session #2


cd('R:\ASommer\PilocarpineRecordings\PTP_5.1\2-11-2026(9.8)\RHS_260211_091012') % change current folder

v = load('auto_sz_classifier.mat');  % load model
v1 = load('PTP 5.1.sessiondata.mat'); % load ground truth; ensure correct animalID/file name

ts  = 1:size(v.estimateLabel,1);    % time vector
% 
sz_ep = [cell2mat({v1.sessiondata.events.start}') cell2mat({v1.sessiondata.events.stop}')]; % load Ground Truth start and end pts for confirmed seizure events

kp = InIntervals(ts,sz_ep);     % kp == GT sz indices
    % tests which values fall in a list of intervals (start/stop pairs) == logical response

auto_sz =  v.estimateLabel(:,1)==2; % events detected and classified by model as sz
    % class 2 = seizure

ons = find(diff([0;auto_sz])>0);    % find detected sz "ons"
offs = find(diff([0;auto_sz])<0);   % find detected sz "offs"


% if recording stopped mid-seizure, offs = end of recording
if length(ons) > length(offs)
    maxT = sm_getFileDur(v1.sessiondata.lfp_file); % from GT data (v1)
    offs = [offs;maxT];
end

% find detected sz events >10s
sz_ep_auto = [ons offs];            
kp_auto = diff(sz_ep_auto,[],2)>10;

sz_ep_auto = sz_ep_auto(kp_auto,:);
kp_auto = InIntervals(ts,sz_ep_auto);

%detected as seizure
seiz_auto_short = v.estimateLabel(:,1)==2 & ~kp_auto; % logical operation == detected sz AND NOT detected sz >10s

scores = v.estimateLabel(:,3); % col 3 =  prob from cat 2 = seizure

scores = scores(:);   % ensure column vectors
scoresDur = scores;   % ensure column vectors


scoresDur(seiz_auto_short) = 0; % hard code scores of short auto seizures to 0

labels = double(kp); % type cast binary into double


% Basic validation
assert(numel(scores) == numel(labels), ...
    'Scores and labels must be the same length.');



%% Test 3 ManualDetect Import


% cd('R:\ASommer\PilocarpineRecordings\PTP_5.1\2-11-2026(9.8)\RHS_260211_091012') % change current folder
cd('R:\ASommer\PilocarpineRecordings\Tests\UnpluggedTest\2-17-2026(14.45)\RHS_260217_145001')
v = load('auto_sz_classifier_refactor.mat');  % load model


ts  = 1:size(v.estimateLabel,1);    % time vector


ev = LoadEvents('manualDetect.evt.szr');
sz_starts=ev.time([contains(ev.description,'on') & contains(ev.description,'sz')]);
sz_ends=ev.time([contains(ev.description,'off') & contains(ev.description,'sz')]);
sz_times=[sz_starts, sz_ends];
sz_truth = InIntervals(ts,sz_times);     % kp == GT sz indices

noise_starts=ev.time([contains(ev.description,'on') & contains(ev.description,'noise')]);
noise_ends=ev.time([contains(ev.description,'off') & contains(ev.description,'noise')]);

if ~isempty(noise_starts) && ~isempty(noise_ends)
    noise_times=[noise_starts, noise_ends];
    noise_truth = 2*InIntervals(ts,noise_times)';     % kp == GT sz indices
else
    noise_truth=zeros(size(ts))';
end

true_labels=double(sz_truth)+double(noise_truth)+1;

class_labels =  v.estimateLabel(:,1);


auto_sz =  v.estimateLabel(:,1)==2; % events detected and classified by model as sz
    % class 2 = seizure

ons = find(diff([0;auto_sz])>0);    % find detected sz "ons"
offs = find(diff([0;auto_sz])<0);   % find detected sz "offs"


% if recording stopped mid-seizure, offs = end of recording
if length(ons) > length(offs)
    maxT = sm_getFileDur(v1.sessiondata.lfp_file); % from GT data (v1)
    offs = [offs;maxT];
end

% find detected sz events >10s
sz_ep_auto = [ons offs];            
kp_auto = diff(sz_ep_auto,[],2)>10;

sz_ep_auto = sz_ep_auto(kp_auto,:);
kp_auto = InIntervals(ts,sz_ep_auto);

%detected as seizure
scores = v.estimateLabel(:,3); % col 3 =  prob from cat 2 = seizure

scores = scores(:);   % ensure column vectors
scoresDur = scores;   % ensure column vectors

seiz_auto_short = v.estimateLabel(:,1)==2 & ~kp_auto; % logical operation == detected sz AND NOT detected sz >10s
scoresDur(seiz_auto_short) = 0; % hard code scores of short auto seizures to 0

% Basic validation
assert(numel(scores) == numel(true_labels), ...
    'Scores and labels must be the same length.');


%% Plot model confidence probabilities over time (? unsure of official plot name)

figure;
plot(ts,kp_auto*20,'Color',"r");
hold on
plot(ts,true_labels*22,'Color',"k");

plot(ts,scores,'Color',"b");
plot(ts,scoresDur,'Color',"m");
%%

cm=confusionmat(true_labels,class_labels);
cm=confusionchart(true_labels,class_labels);

%% Compute and plot ROC curve

% X = false positive rate
% Y = true positive rate
% T = decision thresholds
% AUC = area under the ROC curve

[Xdur, Ydur, Tdur, AUCdur] = perfcurve(true_labels, scoresDur, 1);
for i = unique(true_labels)
    display(i)
    scores = v.estimateLabel(:,i+1); % col 3 =  prob from cat 2 = seizure
    [X, Y, T, AUC] = perfcurve(true_labels==i, scores, 1);
    
    % Plot ROC
    figure;
    plot(X, Y, 'LineWidth', 2);
    hold on;
    plot(Xdur, Ydur, 'LineWidth', 2);
    plot([0 1], [0 1], 'k--');   % chance line
    hold off;
    
    xlabel('False Positive Rate');
    ylabel('True Positive Rate');
    title(sprintf('ROC Curve, Class %1i, AUC = %.3f', i, AUC));
    grid on;
end

%%

% cd('R:\ASommer\PilocarpineRecordings\PTP_5.1\2-11-2026(9.8)\RHS_260211_091012') % change current folder
cd('R:\ASommer\PilocarpineRecordings\Tests\UnpluggedTest\2-17-2026(14.45)\RHS_260217_145001')
v = load('auto_sz_classifier_refactor.mat');  % load model


ts  = 1:size(v.estimateLabel,1);    % time vector


ev = LoadEvents('manualDetect.evt.szr');
sz_starts=ev.time([contains(ev.description,'on') & contains(ev.description,'sz')]);
sz_ends=ev.time([contains(ev.description,'off') & contains(ev.description,'sz')]);
sz_times=[sz_starts, sz_ends];
sz_truth = InIntervals(ts,sz_times);     % kp == GT sz indices

noise_starts=ev.time([contains(ev.description,'on') & contains(ev.description,'noise')]);
noise_ends=ev.time([contains(ev.description,'off') & contains(ev.description,'noise')]);

if ~isempty(noise_starts) && ~isempty(noise_ends)
    noise_times=[noise_starts, noise_ends];
    noise_truth = 2*InIntervals(ts,noise_times)';     % kp == GT sz indices
else
    noise_truth=zeros(size(ts))';
end

true_labels=double(sz_truth)+double(noise_truth)+1;

class_labels =  v.estimateLabel(:,1);

cm=confusionchart(true_labels,class_labels);

roc=rocmetrics(true_labels, v.estimateLabel(:,2:3),[1,2]);
plot(roc)