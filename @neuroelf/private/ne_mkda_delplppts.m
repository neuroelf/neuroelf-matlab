% PUBLIC FUNCTION ne_mkda_delplppts: delete PLP points
function varargout = ne_mkda_delplppts(varargin)

% Version:  v0.9c
% Build:    11120710
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

% get selection
if nargin > 2 && ...
    isa(varargin{3}, 'double') && ...
   ~isempty(varargin{3}) && ...
    numel(varargin{3}) == max(size(varargin{3})) && ...
    ndims(varargin{3}) < 3 && ...
   ~any(isinf(varargin{3}) | isnan(varargin{3}) | varargin{3} < 1 | varargin{3} > size(plp.Points, 1))
    pidx = unique(fix(varargin{3}(:)));
else
    pidx = ch.Points.Value;
    if ~isempty(pidx)
        pidx(pidx < 3) = [];
    end
    if ~isempty(pidx)
        pidx = ch.Points.UserData{2}(pidx(:) - 2);
    end
end

% remove points
if ~isempty(pidx)
    plp.Points(pidx, :) = [];
end

% update listed points
lbt = ch.Points.ListboxTop;
ne_mkda_listpoints;
ch.Points.ListboxTop = min(lbt, length(ch.Points.String));
drawnow;

% store saved flag
plp.RunTimeVars.Saved = false;
