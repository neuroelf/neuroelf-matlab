function [output, status] = ne_freesurfer(command, opts)
% ne_freesurfer  - call FreeSurfer from Matlab
%
% FORMAT:       [output, status] = ne_freesurfer(command [, opts])
%
% Input fields:
%
%       command     one of
%                   'recon-all'
%       opts        optional settings (depending on command)
%        .subdir    SUBJECTS_DIR content (set to pwd if not set)
%
% Output fields:
%
%       output      text output returned by command
%       status      return (exit) status of command
%
% Note: this is merely a wrapper that adds the configured FreeSurferHome
%       from freesurfer.ini to the path and then runs the command

% persistent ini file
persistent fsini;
if numel(fsini) ~= 1 || ...
   ~isstruct(fsini)

    % load ini file
    fsiniobj = xini([neuroelf_path('config') '/freesurfer.ini'], 'convert');
    fsini = fsiniobj.GetComplete;
    fsiniobj.Release;

    % store additional information
    fsini.OS = regexprep(computer, '\d+', '');
    fsini.FSHome = fsini.FreeSurferHome.(fsini.OS);

    % shell
    if isempty(regexpi(fsini.OS, '^win'))
        try
            [shfld, fsini.Shell] = fileparts(invsystem('echo $0'));
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            fsini.Shell = 'bash';
        end
    else
        fsini.Shell = 'C:\Windows\System32\cmd.exe';
    end
end

% argument check
if nargin < 1 || ...
   ~ischar(command) || ...
    isempty(command) || ...
    isempty(regexpi(command(:)', '^s*(recon-all)\s+'))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing command argument.' ...
    );
end
if nargin < 2 || ...
    numel(opts) ~= 1 || ...
   ~isstruct(opts)
    opts = struct;
end
if ~isfield(opts, 'subdir') || ...
   ~ischar(opts.subdir) || ...
    isempty(opts.subdir) || ...
    exist(opts.subdir(:)', 'dir') ~= 7
    opts.subdir = pwd;
else
    opts.subdir = opts.subdir(:)';
end

% run command
try
    switch (fsini.Shell)
        case {'bash'}
            [output, status] = invsystem( ...
                sprintf('export FREESURFER_HOME=%s ; export SUBJECTS_DIR=%s ; source $FREESURFER_HOME/SetUpFreeSurfer.sh ; %s', ...
                fsini.FSHome, opts.subdir, command));
        otherwise
            error( ...
                'neuroelf:UnsupportedShell', ...
                'Unsupported shell: %s.', ...
                fsini.Shell ...
            );
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    output = ne_eo.message;
    status = 255;
end
