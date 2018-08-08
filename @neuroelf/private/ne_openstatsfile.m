% FUNCTION ne_openstatsfile: open file as stats file
function varargout = ne_openstatsfile(varargin)
% ne_openstatsfile  - open a file as a statistical object
%
% FORMAT:       [obj = ] ne_openstatsfile(src, evt, filespec)
%
% Input fields:
%
%       src         handle of UI object issuing the call (0 for GUI)
%       evt         event data (currently unused)
%       filespec    either filename or object reference (optional)
%
% Output fields:
%
%       obj         object reference (e.g. for filenames)

% Version:  v1.1
% Build:    16061711
% Date:     Jun-17 2016, 11:56 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2015, 2016, Jochen Weber
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of Columbia University nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% global variable
global ne_gcfg;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% not while in a callback
if ne_gcfg.c.incb
    return;
end
ne_gcfg.c.incb = true;

% get handles
ch = ne_gcfg.h;

% current status of pointer
mfp = ch.MainFig.Pointer;
ch.MainFig.Pointer = 'watch';

% filename given and valid
if nargin > 2 && ...
    ischar(varargin{3}) && ...
   ~isempty(varargin{3}) && ...
    exist(varargin{3}, 'file') > 0 && ...
   ~isempty(regexpi(varargin{3}(:)', '\.(hdr|head|img|nii)$'))
    hname = varargin{3}(:)';
    if regexpi(hname, '\.img$')
        hname = regexprep(hname, '\.img$', '.hdr', 'preservecase');
    end

% object given and valid
elseif nargin > 2 && ...
    numel(varargin{3}) == 1 && ...
    isxff(varargin{3}, {'hdr', 'head'})
    hname = varargin{3};
    
% otherwise request
else
    [hname, hpath] = uigetfile( ...
        {'*.hdr;*.nii', 'Analyze files (*.hdr, *.nii)'; ...
         '*.nii.gz', 'Compressed Analyze files (*.nii.gz)'; ...
         '*.head', 'AFNI HEAD files (*.head)'}, ...
         'Please select a statistical file to load...', ...
         'MultiSelect', 'on');
    if isequal(hname, 0) || ...
        isequal(hpath, 0) || ...
        isempty(hname)
        ch.MainFig.Pointer = mfp;
        ne_gcfg.c.incb = false;
        return;
    end
    if isempty(hpath)
        hpath = pwd;
    elseif any(hpath(end) == '\/')
        hpath(end) = [];
    end
    hpath = strrep(hpath, filesep, '/');
    if ischar(hname)
        hname = {hname};
    end
    for fc = numel(hname):-1:1
        hname{fc} = [hpath '/' hname{fc}];
        if exist(hname{fc}, 'file') ~= 2
            hname(fc) = [];
        end
    end
    if isempty(hname)
        ch.MainFig.Pointer = mfp;
        ne_gcfg.c.incb = false;
        return;
    end
end

% load file
try
    if numel(hname) == 1
        if iscell(hname)
            hobj = xff(hname{1});
        else
            hobj = hname;
        end
        if ~isxff(hobj, {'hdr', 'head'})
            clearxffobjects({hobj});
            error( ...
                'neuroelf:BadFileSelected', ...
                'Invalid file selected.' ...
            );
        end
    else
        hobj = xff(hname);
        if numel(hobj) == 1
            if ~isxff(hobj, {'hdr', 'head'})
                clearxffobjects({hobj});
                error( ...
                    'neuroelf:BadFileSelected', ...
                    'Invalid file selected.' ...
                );
            end
        else
            clearxffobjects(hobj);
            error( ...
                'neuroelf:BadFileSelected', ...
                'Invalid file selected.' ...
            );
        end
    end
catch ne_eo;
    ch.MainFig.Pointer = mfp;
    ne_gcfg.c.incb = false;
    rethrow(ne_eo);
end
if nargout > 0
    varargout{1} = hobj;
end

% ensure MapSelection
if ~isfield(hobj.RunTimeVars, 'MapSelection')
    hobj.RunTimeVars.MapSelection = {{}, []};
end

% set stats object status
hobj.RunTimeVars.StatsObject = true;

% request missing information (distribution type, DF)
if nargin < 4 || ...
   ~iscell(varargin{4}) || ...
    numel(varargin{4}) ~= 2 || ...
   ~ischar(varargin{4}{1}) || ...
    numel(varargin{4}{1}) ~= 1 || ...
   ~any(lower(varargin{4}{1}) == 'ftr%') || ...
   ~ischar(varargin{4}{2}) || ...
    isempty(regexpi(varargin{4}{2}(:)', '^\d+$'))
    ra = inputdlg({'Stats type (t, F, r, %):', 'Degrees of freedom:'}, ...
        'NeuroElf - user input', 1, {'  t', '  20'});
else
    ra = varargin{4};
end
if ~iscell(ra) || ...
    numel(ra) ~= 2 || ...
   ~ischar(ra{1}) || ...
   ~ischar(ra{2}) || ...
    isempty(ra{1}) || ...
    isempty(ra{2})
    ch.MainFig.Pointer = mfp;
    clearxffobjects({hobj});
    ne_gcfg.c.incb = false;
    return;
end
ra{1} = lower(ddeblank(ra{1}));
if numel(ra{1}) ~= 1 || ...
   ~any('ftr%' == ra{1})
    ch.MainFig.Pointer = mfp;
    clearxffobjects({hobj});
    ne_gcfg.c.incb = false;
    return;
end
ra{2} = ddeblank(ra{2});
try
    ra{2} = eval(['[' ra{2} ']']);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    ra{2} = [234, 1];
end

% set map status
rtv = hobj.RunTimeVars;
mnames = {rtv.Map.Name};
if isempty(mnames) || ~all(strcmpi(mnames{1}, mnames(:)))
    hname = {};
end    
for mc = 1:numel(rtv.Map)
    rtv.Map(mc).DF1 = ra{2}(1);
    if numel(ra{2}) > 1
        rtv.Map(mc).DF2 = ra{2}(2);
    else
        rtv.Map(mc).DF2 = 1;
    end
    switch(ra{1})
        case {'f'}
            rtv.Map(mc).Type = 4;
        case {'r'}
            rtv.Map(mc).Type = 2;
            rtv.Map(mc).LowerThreshold = 0.1;
            rtv.Map(mc).UpperThreshold = 0.5;
        case {'t'}
            rtv.Map(mc).Type = 1;
        case {'%'}
            rtv.Map(mc).Type = 15;
    end
    if iscell(hname) && numel(hname) == numel(rtv.Map)
        [hpath, hfname] = fileparts(hname{mc}(:)');
        if ~isempty(hfname)
            rtv.Map(mc).Name = hfname;
        end
    end
end
hobj.RunTimeVars = rtv;

% open in NeuroElf
ch.MainFig.Pointer = mfp;
ne_gcfg.c.incb = false;
ne_openfile(0, 0, hobj, true);
