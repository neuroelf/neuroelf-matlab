function createprtsfromeprime(efile, opts)
% createprtsfromeprime  - create PRT (Protocol) files from EPrime text file
%
% FORMAT:       createprtsfromeprime([efile [, opts]])
%
% Input fields:
%
%       efile       filename of EPrime output file
%       opts        optional settings
%        .cond      1xC sub-struct with condition rules
%         .color    1x3 RGB code
%        .runfield  fieldname from which to detect runs
%        .syncfield fieldname from which to get scanner sync
%
% No output fields.

% for now use all of neuroelf
using(neuroelf, 'all');

% argument check
if nargin < 1 || ...
   ~ischar(efile) || ...
    exist(efile(:)', 'file') ~= 2
    [efile, epath] = uigetfile( ...
        {'*.txt', 'EPrime log file exports (*.txt)'; ...
         '*.xls', 'EPrime log file Excel sheet (*.xls)'; ...
         '*.*'  , 'All files (*.*)'}, ...
        'Please select an EPrime log file...');
    if isequal(efile, 0) || ...
        isequal(epath, 0)
        return;
    end
    if isempty(epath)
        epath = pwd;
    end
    efile = fullfile(epath, efile);
end
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'syncfield') || ...
   ~ischar(opts.syncfield) || ...
    isempty(opts.syncfield) || ...
   ~isvarname(opts.syncfield(:)')
    opts.syncfield = '';
else
    opts.syncfield = opts.syncfield(:)';
end
if ~isfield(opts, 'runfield') || ...
   ~ischar(opts.runfield) || ...
    isempty(opts.runfield) || ...
   ~isvarname(opts.runfield(:)')
    opts.runfield = opts.syncfield;
else
    opts.runfield = opts.runfield(:)';
end

% try to read EPrime file
try
    edat = readeprimetextlog(efile);
    log = edat.Log;
    % handle multiple levels
    if numel(log) < 4
        error( ...
            'neuroelf:BadFileContent', ...
            'Too few events to create protocol.' ...
        );
    end
catch ne_eo;
    rethrow(ne_eo);
end

% check content
if ~isfield(log, opts.syncfield) || ...
   ~isfield(log, opts.runfield)
    error( ...
        'neuroelf:BadArgument', ...
        'Log must contain sync and run fields.' ...
    );
end
runs = cat(1, log.(opts.runfield));
urun = unique(runs, 'stable');
nrun = numel(urun);

% iterate over runs

