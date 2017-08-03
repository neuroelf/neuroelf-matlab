function [varargout] = mtc_concat(targetfile, mtcs, varargin)
% mtc_concat  - concatenate MTCs
%
% FORMAT:       newmtc = mtc_concat(targetfile, mtclist [, options]);
%         or    newmtc = mtc_concat(targetfile, mtc1, mtc2, ... [, options]);
%
% Input fields:
%
%       targetfile  filename of MTC to write
%       mtclist     cell array with MTCs to concatenated
%       mtc1, ...   single MTCs to concatenated
%       options     1x1 struct with optional settings
%        .trans     transformation to apply on each MTC on reading
%                   either of 'psc' or 'z'
%
% Output fields:
%
%       newmtc      object handle to newly written MTC

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:28 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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

% argument check
if nargin < 2 || ...
   ~ischar(targetfile) || ...
    isempty(targetfile) || ...
   (~all(isxff(mtcs, 'mtc')) && ...
    ~iscell(mtcs)) || ...
    isempty(mtcs)
    error( ...
        'neuroelf:BadArgument', ...
        'You must give a valid target filename and some MTCs.' ...
    );
end
if all(isxff(mtcs, 'mtc'))
    mtcc = cell(1, numel(mtcs));
    for mc = 1:numel(mtcs)
        mtcc{mc} = mtcs(mc);
    end
    mtcs = mtcc;
end
for n = 3:nargin
    if isxff(varargin{n - 2}, 'mtc')
        mtcs{end + 1} = varargin{n - 2};
    end
end
nm = numel(mtcs);
if nm == 1
    error( ...
        'neuroelf:TooFewArguments', ...
        'You need at least two MTCs to concatenate.' ...
    );
end

% options ?
transopt = 0;
if nargin > 2 && ...
    isstruct(varargin{end}) && ...
    numel(varargin{end}) == 1
    options = varargin{end};
else
    options = struct;
end
if isfield(options, 'trans') && ...
    ischar(options.trans) && ...
   ~isempty(options.trans)
    switch (lower(options.trans(:)'))
        case {'psc', '%'}
            transopt = 1;
        case {'z'}
            transopt = 2;
    end
end

% check mtcs argument for filenames
loadmtcs = false(1, nm);
for mc = 1:nm
    if ischar(mtcs{mc}) && ...
        exist(mtcs{mc}(:)', 'file') == 2
        loadmtcs(mc) = true;
    elseif numel(mtcs{mc}) ~= 1 || ...
       ~isxff(mtcs{mc}, 'mtc')
        error( ...
            'neuroelf:BadArgument', ...
            'Cell array must contain either MTC objects or filenames.' ...
        );
    end
end
xffroot = xff();
if any(loadmtcs)
    mtctiosz = xffroot.TransIOSize('mtc', 65536);
    rlsame = xffroot.Config('reloadsame', true);
    for mc = find(loadmtcs)
        if ischar(mtcs{mc})
            try
                mtcfname = mtcs{mc}(:)';
                mtcs{mc} = xff(mtcfname);
            catch ne_eo;
                clearxffobjects(mtcs(loadmtcs));
                xffroot.Config('reloadsame', rlsame);
                xffroot.TransIOSize('mtc', mtctiosz);
                error( ...
                    'neuroelf:xffError', ...
                    'Error opening MTC ''%s'' (%s).', ...
                    mtcfname, ne_eo.message ...
                );
            end
        end
    end
    xffroot.Config('reloadsame', rlsame);
    xffroot.TransIOSize('mtc', mtctiosz);
end

% check headers
mh = mtcs{1};
vs = size(mh.MTCData);
mtp = zeros(1, nm);
ntp = vs(1);
mtp(1) = ntp;
for mc = 2:nm
    cs = size(mtcs{mc}.MTCData);
    if cs(2) ~= vs(2)
        clearxffobjects(mtcs(loadmtcs));
        error( ...
            'neuroelf:BadArgument', ...
            'MTCs must match in NrOfVertices property.' ...
        );
    end
    ntp = ntp + cs(1);
    mtp(mc) = cs(1);
end
mtp = 1 + [0, cumsum(mtp)];

% construct transioobject
tobjname = tempname;
tobj = transio(tobjname, 'ieee-le', 'single', 0, [ntp, cs(2)], 1);

% copy first MTC object, and set new NrOfTimePoints, MTCData
newmtc = mtcs{1}.CopyObject;
newmtc.NrOfTimePoints = ntp;
newmtc.SourceVTCFile = '<concatenated>';
newmtc.LinkedPRTFile = '<none>';
newmtc.MTCData = tobj;

% create temporary object in mem
mo = single(0);
mo(ntp, 1000) = 0;

% always take 1000 vertices at a time
sc = 1;
st = 999;
while (sc + st) <= cs(2)
    for mc = 1:nm
        if transopt == 2
            mo(mtp(mc):(mtp(mc + 1)-1), :) = ...
                ztrans(double(mtcs{mc}.MTCData(:, sc:(sc + st))));
        elseif transopt == 1
            mo(mtp(mc):(mtp(mc + 1)-1), :) = ...
                psctrans(double(mtcs{mc}.MTCData(:, sc:(sc + st))));
        else
            mo(mtp(mc):(mtp(mc + 1)-1), :) = ...
                mtcs{mc}.MTCData(:, sc:(sc + st));
        end
    end
    newmtc.MTCData(:, sc:(sc + st)) = mo;
    sc = sc + st + 1;
end

% for the last part, remove part of mo
mo(:, (cs(2) + 2 - sc):end) = [];
if ~isempty(mo)
    for mc = 1:nm
        if transopt == 2
            mo(mtp(mc):(mtp(mc + 1)-1), :) = ...
                ztrans(double(mtcs{mc}.MTCData(:, sc:cs(2))));
        elseif transopt == 1
            mo(mtp(mc):(mtp(mc + 1)-1), :) = ...
                psctrans(double(mtcs{mc}.MTCData(:, sc:cs(2))));
        else
            mo(mtp(mc):(mtp(mc + 1)-1), :) = ...
                mtcs{mc}.MTCData(:, sc:cs(2));
        end
    end
    newmtc.MTCData(:, sc:cs(2)) = mo;
end

% clear objects
clearxffobjects(mtcs(loadmtcs));

% save MTC under new name
try
    newmtc.SaveAs(targetfile);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    warning( ...
        'neuroelf:BadArgument', ...
        'Bad target filename given. Please remove ''%s'' manually.', ...
        tobjname ...
    );
    newmtc.ClearObject;
    return;
end

% delete old temp object
delete(tobjname);

% if no argout, clear object
if nargout < 1
    newmtc.ClearObject;
else
    varargout = cell(1, nargout);
    varargout{1} = newmtc;
end
