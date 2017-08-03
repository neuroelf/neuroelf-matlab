% use neuroelf
n = neuroelf;

% select file
s = uigetfile('*.xls');
if isequal(s, 0) || ...
    isempty(s)
    error('No file selected.');
end

% generate onsets
onsets = struct;

% generate output filename
o = strrep(s, '.xls', '_onsets.mat');

% load file
e = n.readeprimetextlog(s);

% load subjects list
load slist;

% get subjects, and create unique list
subjects = cat(1, e.Log.Subject);
usubs = unique(subjects);

% iterate over subjects
for sc = 1:numel(usubs)
    
    % generate data matrix
    dmtx = NaN .* zeros(108, 5);
    
    % find rows for subject
    subrows = find(subjects == usubs(sc));
    
    % subject ID
    if any(slist(:, 2) == usubs(sc))
        sid = sprintf('s%04d', slist(slist(:, 2) == usubs(sc), 1));
    else
        sid = sprintf('s1%03d', usubs(sc));
    end
    
    % iterate over sessions
    for ssc = 1:2
        
        % find rows matching this
        sessrows = subrows(cat(1, e.Log(subrows).Session) == ssc);
        
        % no data, continue
        if isempty(sessrows)
            continue;
        end
        
        % get data
        sessdata = e.Log(sessrows);
        
        % iterate over runs
        for rc = 1:3
            
            % not enough data
            if numel(sessdata) < (18 * rc)
                break;
            end
            
            % get rundata
            rundata = sessdata(18*rc-17:18*rc);
            
            % target rows
            trows = 54 * (ssc - 1) + 18 * (rc - 1) + 1;
            trows = trows:(trows+17);
            
            % get run onset
            firststimdelay = rundata(1).StimSlide_OnsetDelay;
            runonset = rundata(1).StimSlide_OnsetTime - (23000 + firststimdelay);
            
            % get stim condition
            stimcond = {rundata.conditionname};
            stimcond(~cellfun('isempty', regexpi(stimcond(:), 'neu'))) = {'1'};
            stimcond(~cellfun('isempty', regexpi(stimcond(:), 'look.*neg'))) = {'2'};
            stimcond(~cellfun('isempty', regexpi(stimcond(:), 'rea'))) = {'3'};
            dmtx(trows, 1) = str2double(stimcond(:));
            
            % get cue, stim, and rating onsets, and rating values
            dmtx(trows, 2) = cat(1, rundata.StimSlide_OnsetTime) - (runonset + 2000 + firststimdelay);
            dmtx(trows, 3) = cat(1, rundata.StimSlide_OnsetTime) - runonset;
            dmtx(trows, 4) = cat(1, rundata.RatingSlide_OnsetTime) - runonset;
            dmtx(trows, 5) = cat(1, rundata.RatingSlide_RESP);
        end
    end
    
    % store matrix
    onsets.(sid) = dmtx;
end
