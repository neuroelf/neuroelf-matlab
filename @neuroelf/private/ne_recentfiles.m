% FUNCTION ne_recentfiles: update list of recently opened files
function ne_recentfiles(varargin)

% Version:  v0.9b
% Build:    11101810
% Date:     Aug-11 2010, 9:00 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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

% only once at a time
if ne_gcfg.c.incb
    return;
end
ne_gcfg.c.incb = true;

ch = ne_gcfg.h.RecentFiles;
rf = ne_gcfg.c.ini.RecentFiles;
rff = varargin{3};
rif = varargin{4};

% delete unavailable
if rf.DeleteUnavailable > 0
    for fc = numel(rf.(rif)):-1:1
        if exist(rf.(rif){fc}, 'file') ~= 2
            rf.(rif)(fc) = [];
        end
    end
end

% set initial visible stat
for fc = 1:rf.Number
    set(ch.(rff){fc}, 'Visible', 'off');
end

% function handles
fh_openfile = @ne_openfile;
fh_openstatsfile = @ne_openstatsfile;

% add files
if strcmp(rif, 'StatsVar')
    for fc = 1:numel(rf.(rif))
        fname = rf.(rif){fc};
        if numel(fname) > 99
            flabel = [fname(1:47) ' ... ' fname(end-46:end)];
        else
            flabel = fname;
        end
        if isempty(regexpi(fname, '\.(hdr|nii)$'))
            set(ch.(rff){fc}, ...
                'Callback', {fh_openfile, fname}, ...
                'Label',    flabel, ...
                'Visible',  'on');
        else
            set(ch.(rff){fc}, ...
                'Callback', {fh_openstatsfile, fname}, ...
                'Label',    flabel, ...
                'Visible',  'on');
        end
    end
else
    for fc = 1:numel(rf.(rif))
        fname = rf.(rif){fc};
        if numel(fname) > 99
            flabel = [fname(1:47) ' ... ' fname(end-46:end)];
        else
            flabel = fname;
        end
        set(ch.(rff){fc}, ...
            'Callback', {fh_openfile, fname}, ...
            'Label',    flabel, ...
            'Visible',  'on');
    end
end

% remove only-once flag
ne_gcfg.c.incb = false;
