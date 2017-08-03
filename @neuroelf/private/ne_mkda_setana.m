% PUBLIC FUNCTION ne_mkda_setana: set UIs to current analysis' values
function varargout = ne_mkda_setana(varargin)

% Version:  v0.9c
% Build:    13012611
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
hFig = ne_gcfg.h.MKDA.MKDAFig;
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
rtv = plp.RunTimeVars;
anas = rtv.MKDAAnalyses;

% inputs
if nargin < 3 || ...
   ~isa(varargin{3}, 'double') || ...
    numel(varargin{3}) ~= 1 || ...
    isinf(varargin{3}) || ...
    isnan(varargin{3}) || ...
    varargin{3} < 1 || ...
    varargin{3} > size(anas, 1) || ...
    varargin{3} ~= fix(varargin{3})
    anaidx = ch.Analyses.Value;
else
    anaidx = varargin{3};
end
ch.Analyses.Value = anaidx;

% update condition, contrast, weighting, and statistical unit
ana = anas{anaidx, 2};
cntcol = findfirst(strcmpi(ana.ContColumn, ch.ContColumn.String));
if isempty(cntcol)
    ch.ContColumn.Value = 1;
    ne_mkda_setcontcol;
else
    ch.Contrast.String = ana.Contrast;
    ch.ContColumn.Value = cntcol;
end
ch.CndParts.Value = [];
ch.CndParts.String = ana.CndParts;
ch.Weights.String = ana.Weights;
stdcol = findfirst(strcmpi(ana.StudyColumn, ch.StudyColumn.String));
if ~isempty(stdcol)
    ch.StudyColumn.Value = stdcol;
end

% update general settings
ch.Iterations.String = sprintf('%d', ana.Iterations);
mskstring = ch.Mask.String;
if ~iscell(mskstring)
    mskstring = cellstr(mskstring);
end
mskstring{1} = ana.Mask;
ch.Mask.String = mskstring;
ch.Mask.Value = 1;
ch.Resolution.Value = ana.Resolution;
if ana.Scaling(1) == 'i'
    hFig.RadioGroupSetOne('Scale', 1);
    hFig.SetGroupEnabled('Gauss', 'off');
else
    hFig.RadioGroupSetOne('Scale', 2);
    hFig.SetGroupEnabled('Gauss', 'on');
end
ch.SphereSize.String = sprintf('%.4g', ana.SphereSize);
ch.SphereTaper.String = sprintf('%.4g', ana.SphereTaper);
if ana.NullDist(1) == 's'
    hFig.RadioGroupSetOne('SUNull', 1);
else
    hFig.RadioGroupSetOne('SUNull', 2);
end

% parse contrast
ne_mkda_parsecont;

% list points
ne_mkda_listpoints;
