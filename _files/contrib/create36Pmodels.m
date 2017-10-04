% script that takes an MDM file (file selector) with PRTs
% and creates SDMs that contain configurable regressors

% use NeuroElf
n = neuroelf;

% load MDM
mdm = xff('*.mdm', 'Please select MDM file with PRT references...');

% check MDM
if ~isxff(mdm, 'mdm') || isempty(mdm.XTC_RTC) || ...
   ~isfield(mdm.RunTimeVars, 'MotionParameters') || ...
    numel(mdm.RunTimeVars.MotionParameters) ~= size(mdm.XTC_RTC, 1) || ...
    any(cellfun('isempty', regexpi(mdm.XTC_RTC(:, 2), '\.prt$')))
    if isxff(mdm)
        mdm.ClearObject;
    end
    error('neuroelf:scriptError', 'Invalid MDM file.');
end

% get options
[nuisopt, nselok] = listdlg( ...
    'PromptString', 'Please select the nuisance regressors added to the SDMs:', ...
    'SelectionMode', 'multiple', ...
    'ListString', {'MotionParameters'; 'MotionParameters derivative'; ...
        'MotionParameters squared'; 'MotionParameters squared derivative'; ...
        'GlobalSignal'; 'GlobalSignal derivative'; ...
        'GlobalSignal squared'; 'GlobalSignal squared derivative'; ...
        'SpikeRegressors'; 'TemporalFilters'}, ...
    'InitialValue', (1:9)');
if ~isequal(nselok, 1) || isempty(nuisopt)
    mdm.ClearObject;
    error('neuroelf:scriptError', 'Nothing selected.');
end

% options
motpars = any(nuisopt == 1);
motparsd = any(nuisopt == 2);
motparsq = any(nuisopt == 3);
motparsqd = any(nuisopt == 4);
globsig = any(nuisopt == 5);
globsigd = any(nuisopt == 6);
globsigq = any(nuisopt == 7);
globsigqd = any(nuisopt == 8);
spkregs = any(nuisopt == 9);
tempfilt = any(nuisopt == 10);

% ask for global signal masks
% globsigs = uigetfile(...);

% load global signal files

% get PRT files
vtcs = mdm.XTC_RTC(:, 1);
prts = mdm.XTC_RTC(:, 2);

% iterate over data
for rc = 1:numel(prts)
    
    % load VTC (transio) and PRT
    vtc = xff(vtcs{rc}, 't');
    prt = xff(prts{rc});
    if ~isxff(vtc, 'vtc') || ~isxff(prt, 'prt')
        if isxff(vtc)
            vtc.ClearObject;
        end
        if isxff(prt)
            prt.ClearObject;
        end
        error('neuroelf:scriptError', 'Files for run %d not VTC/PRT.', rc);
    end

    % 
end
