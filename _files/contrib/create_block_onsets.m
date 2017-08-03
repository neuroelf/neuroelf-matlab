function create_block_onsets(filename)

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
matfilename = [filename(1:end-3) 'mat'];

% get different task types (per trial)
if isfield(log.Log, 'Type')
    tasktype = {log.Log.Type};
else
    tasktype = {log.Log.StimType};
end
tasktype = tasktype(:);

% get sync offset
notask = cellfun('isempty', tasktype);
syncoffset = log.Log(findfirst(~notask)).WaitforSigna_FinishTime;

% get onsets, RESP, and RTs (if used as a parameter)
if isfield(log.Log, 'Target_OnsetTime')
    onsettimes = cat(1, log.Log.Target_OnsetTime) - syncoffset;
    targetdurations = cat(1, log.Log.Target_Duration);
else
    onsettimes = cat(1, log.Log.stim_OnsetTime) - syncoffset;
    targetdurations = cat(1, log.Log.stim_Duration);
end
%rts = cat(1, log.Log.Target_RT);
%resp = {log.Log.Target_RESP};
%resp(cellfun('isempty', resp)) = {NaN};

% remove empty entries
tasktype(notask) = [];

% unique tasks
names = unique(tasktype);

% create cell arrays as needed
onsets = cell(1, numel(names));
durations = onsets;

% for each of the conditions
for cc = 1:numel(onsets)

    % find onsets
    onsidx = find(strcmpi(tasktype, names{cc}));

    % find first occurrence in condition
    firstons = onsidx;
    firstons(1 + find(diff(onsidx(:)) == 1)) = [];

    % find last occurrence in condition
    lastons = onsidx;
    lastons(find(diff(onsidx(:)) == 1)) = [];

    % store in onsets and durations
    onsets{cc} = 1e-3 .* onsettimes(firstons);
    durations{cc} = (1e-3 .* (onsettimes(lastons) + targetdurations(lastons))) - onsets{cc};
end

% deblank condition names
names = ddeblank(names(:))';

% save
save(matfilename, 'names', 'onsets', 'durations');
