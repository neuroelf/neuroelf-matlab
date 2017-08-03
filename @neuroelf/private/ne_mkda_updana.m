% PUBLIC FUNCTION ne_mkda_updana: update analysis with current UI values
function varargout = ne_mkda_updana(varargin)

% Version:  v0.9c
% Build:    13012610
% Date:     Nov-08 2011, 1:08 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Jochen Weber
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
ch = ne_gcfg.h.MKDA.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% get content of PLP dropdown
plps = ch.PLPs;
plpud = plps.UserData;
plpid = plps.Value;
try
    plp = plpud{plpid, 3};
    if numel(plp) ~= 1 || ...
       ~isxff(plp, 'plp')
        error( ...
            'neuroelf:GUI:BadPLP', ...
            'Bad PLP object.' ...
        );
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    return;
end
anaidx = ch.Analyses.Value;
ana = plp.RunTimeVars.MKDAAnalyses{anaidx, 2};

% update condition, contrast, weighting, and statistical unit
ana.ContColumn = ch.ContColumn.String{ch.ContColumn.Value};
ana.Contrast = ch.Contrast.String;
ana.CndParts = ch.CndParts.String;
if ~iscell(ana.CndParts)
    if isempty(ana.CndParts)
        ana.CndParts = {};
    else
        ana.CndParts = cellstr(ana.CndParts);
    end
end
ana.Weights = ch.Weights.String;
ana.StudyColumn = ch.StudyColumn.String{ch.StudyColumn.Value};

% update general settings
try
    ana.Iterations = max(1, min(1000000, ceil(str2double(ch.Iterations.String))));
    if ~isa(ana.Iterations, 'double') || ...
        numel(ana.Iterations) ~= 1 || ...
        isinf(ana.Iterations) || ...
        isnan(ana.Iterations)
        ana.Iterations = 5000;
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    ana.Iterations = 5000;
end
ch.Iterations.String = sprintf('%d', ana.Iterations);
mskstring = ch.Mask.String;
if ~iscell(mskstring)
    mskstring = cellstr(mskstring);
end
ana.Mask = mskstring{1};
ana.Resolution = ch.Resolution.Value;
if ch.SphereIndicator.Value > 0
    ana.Scaling = 'indic';
else
    ana.Scaling = 'toone';
end
try
    ana.SphereSize = max(ana.Resolution, min(30, str2double(ch.SphereSize.String)));
    if ~isa(ana.SphereSize, 'double') || ...
        numel(ana.SphereSize) ~= 1 || ...
        isinf(ana.SphereSize) || ...
        isnan(ana.SphereSize)
        ana.SphereSize = 12;
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    ana.SphereSize = 12;
end
ch.SphereSize = sprintf('%.4g', ana.SphereSize);
try
    ana.SphereTaper = max(0, min(ana.SphereSize, str2double(ch.SphereTaper.String)));
    if ~isa(ana.SphereTaper, 'double') || ...
        numel(ana.SphereTaper) ~= 1 || ...
        isinf(ana.SphereTaper) || ...
        isnan(ana.SphereTaper)
        ana.SphereTaper = 0;
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    ana.SphereTaper = 0;
end
ch.SphereTaper = sprintf('%.4g', ana.SphereTaper);
if ch.NullSpatial.Value > 0
    ana.NullDist = 'spatial';
else
    ana.NullDist = 'units';
end

% update in RunTimeVars
plp.RunTimeVars.MKDAAnalyses{anaidx, 2} = ana;

% save configuration
if ~isempty(plp.FilenameOnDisk)
    try
        plp.SaveRunTimeVars;
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
