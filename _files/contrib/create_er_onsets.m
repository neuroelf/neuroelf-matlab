function create_er_onsets(filename)

% request filename if not given
if nargin < 1 || ...
   ~ischar(filename) || ...
    isempty(filename)
    [filename, filepath] = uigetfile( ...
        {'*.txt', 'Text files (*.txt)'; ...
         '*.log', 'Log files (*.log)'}, 'Please select a log file...');
    if isequal(filename, 0) || ...
        isequal(filepath, 0) || ...
        isempty(filename)
        return;
    end
    if isempty(filepath)
        filepath = pwd;
    end
    filename = [filepath '/' filename];
end

% try to read log file
try
    log = readeprimetextlog(filename, true);
catch ne_eo;
    rethrow(ne_eo);
end

% target file
matfilename = [filename(1:end-4) '.mat'];
matfilename_co = [filename(1:end-4) '_co.mat'];

% get fields that are possible candidates for the sync pulse
logfields = fieldnames(log.Log(1));
numfields = numel(logfields);

% throw out entries that are definitely not helpful
log.Log(cellfun('isempty', {log.Log.Target_ACC})) = [];

% get contents from first log entry
logfieldc = struct2cell(log.Log(1));

% which of these are numbers
logfieldn = cellfun(@isnumeric, logfieldc);

% remove further fields if they change between first and second entry
for fc = 1:numfields
    if ~logfieldn(fc)
        continue;
    end
    if ~isequal(log.Log(1).(logfields{fc}), log.Log(2).(logfields{fc}))
        logfieldn(fc) = false;
    end
end

% which one is the most likely candidate
candidates = cat(2, logfieldc{logfieldn});
bestmatch = candidates(candidates >= 10000);
candidate = findfirst(cellfun(@isequal, logfieldc, repmat({min(bestmatch)}, size(logfieldc))));
if isempty(candidate)
    candidate = findfirst(logfieldn);
end

% map used fields to list
listmap = zeros(numfields, 1);
listmap(logfieldn) = 1:sum(logfieldn);
candidate = listmap(candidate);

% present fields in a list for selection
candidates = logfields(logfieldn);
[syncfield, selok] = listdlg('ListString', candidates, ...
    'SelectionMode', 'single', 'InitialValue', candidate, ...
    'Name', 'Please select the fMRI sync pulse field name...');
if ~isequal(selok, 1)
    return;
end

% get sync offset
syncoffset = log.Log(1).(candidates{syncfield});

% get different pieces of information that are needed
conditiontrigger = cat(1, log.Log.ConditionTrigger);
anticipationonset = cat(1, log.Log.Anticipation_StartTime) - syncoffset;
anticipationduration = cat(1, log.Log.Anticipation_Duration);
feedbackonset = cat(1, log.Log.Feedback_OnsetTime) - syncoffset;
feedbackduration = cat(1, log.Log.Feedback_Duration);
accuracy = cat(1, log.Log.Target_ACC);

% to select correct only, we then add this to the selection
co = (accuracy == 1);

% get the indices for WIN0, WIN1, and WIN2
win0 = (conditiontrigger == 4);
win1 = (conditiontrigger == 3);
win2 = (conditiontrigger == 2);

% create cell arrays as needed
names = {'ant_0', 'ant_1', 'ant_2', 'feed_0', 'feed_1', 'feed_2'};
names_co = {'ant_0_co', 'ant_1_co', 'ant_2_co', 'feed_0_co', 'feed_1_co', 'feed_2_co'};
onsets = cell(1, numel(names));
durations = onsets;

% fill in onsets and durations for full list
onsets{1} = anticipationonset(win0);
onsets{2} = anticipationonset(win1);
onsets{3} = anticipationonset(win2);
onsets{4} = feedbackonset(win0);
onsets{5} = feedbackonset(win1);
onsets{6} = feedbackonset(win2);
durations{1} = anticipationduration(win0);
durations{2} = anticipationduration(win1);
durations{3} = anticipationduration(win2);
durations{4} = feedbackduration(win0);
durations{5} = feedbackduration(win1);
durations{6} = feedbackduration(win2);

% save
save(matfilename, 'names', 'onsets', 'durations');

% then reselect with correct only
onsets{1} = anticipationonset(win0 & co);
onsets{2} = anticipationonset(win1 & co);
onsets{3} = anticipationonset(win2 & co);
onsets{4} = feedbackonset(win0 & co);
onsets{5} = feedbackonset(win1 & co);
onsets{6} = feedbackonset(win2 & co);
durations{1} = anticipationduration(win0 & co);
durations{2} = anticipationduration(win1 & co);
durations{3} = anticipationduration(win2 & co);
durations{4} = feedbackduration(win0 & co);
durations{5} = feedbackduration(win1 & co);
durations{6} = feedbackduration(win2 & co);

% save
names = names_co;
save(matfilename_co, 'names', 'onsets', 'durations');
