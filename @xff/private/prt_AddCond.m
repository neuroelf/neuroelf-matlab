function xo = prt_AddCond(xo, cname, onoffsets, ccolor, pw, pn)
% PRT::AddCond  - add a condition to a PRT file
%
% FORMAT:       [prt = ] prt.AddCond(name, onoffsets [, ccolor [, pw, pn]]);
%
% Input fields:
%
%       name        condition name
%       onoffsets   Ox2 or Ox3 list of on and offsets (and weights)
%       ccolor      condition color
%       pw          parametric weights
%       pn          parametric weight names
%
% Output fields:
%
%       prt         altered PRT object
%
% Note: all input arguments are optional. if none is given, an
%       empty condition with a 'Condition Nr.' identifier and
%       random color is generated
%
% Using: findfirst, varc.

% Version:  v1.1
% Build:    16040710
% Date:     Apr-07 2016, 10:33 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2015, 2016, Jochen Weber
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

% neuroelf library
global ne_methods;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'prt') || ...
   ~ischar(cname) || isempty(cname)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
cname = cname(:)';
pweights = [];
if nargin > 2 && isa(onoffsets, 'double') && ndims(onoffsets) == 2 && size(onoffsets, 2) >= 2 && ...
    all(all(onoffsets(:, 1:2) >= 0)) && ~any(isinf(onoffsets(:)) | isnan(onoffsets(:)))
    if size(onoffsets, 2) > 2
        pweights = onoffsets(:, 3:end);
    end
    onoffsets = round(onoffsets(:, 1:2));
else
    onoffsets = zeros(0, 2);
end
noo = size(onoffsets, 1);
if nargin > 3 && isa(ccolor, 'double') && numel(ccolor) == 3 && ...
   ~any(isinf(ccolor) | isnan(ccolor) | ccolor < 0 | ccolor >= 256)
    ccolor = floor(ccolor(:)');
else
    ccolor = floor(255.99 * rand(1,3));
end
if nargin > 4 && isa(pw, 'double') && ndims(pw) <= 2 && ...
    size(pw, 1) == noo && ~any(isinf(pw(:)))
    pweights = pw;
end
lastp = 0;
if isempty(pweights)
    pweights = zeros(noo, 0);
elseif size(pweights, 1) > 1
    pwvar = (ne_methods.varc(pweights, 1, true) >= 1.5e-8);
    lastp = ne_methods.findfirst(pwvar, -1);
    if isempty(lastp)
        if size(pweights, 2) > 1
            lastp = 0;
        else
            lastp = 1;
        end
    end
    pweights(:, (lastp+1):end) = [];
    if ~isempty(pweights)
        bc.FileVersion = 3;
        bc.ParametricWeights = max(bc.ParametricWeights, size(pweights, 2));
    end
end

% build new condition
newc = struct('ConditionName', {{cname}}, 'NrOfOnOffsets', noo, ...
    'OnOffsets', onoffsets, 'Weights', pweights, 'Color', ccolor);

% parametric weight names
if lastp > 0 && nargin > 5 && iscell(pn) && numel(pn) >= lastp && ...
    all(cellfun(@ischar, pn(pwvar)))
    pn = pn(:)';

    % add to names in RunTimeVars
    if ~isfield(bc.RunTimeVars, 'ParameterNames') || ~iscell(bc.RunTimeVars.ParameterNames)
        bc.RunTimeVars.AutoSave = true;
        bc.RunTimeVars.ParameterNames = cell(1, lastp);
        for pcc = 1:lastp
            bc.RunTimeVars.ParameterNames{1, pcc} = sprintf('p%d', pcc);
        end
    else
        bc.RunTimeVars.ParameterNames = bc.RunTimeVars.ParameterNames(:)';
    end
    if numel(bc.RunTimeVars.ParameterNames) < lastp
        for pcc = 1:lastp
            if pwvar(pcc)
                bc.RunTimeVars.ParameterNames{1, pcc} = sprintf('p%d', pcc);
            end
        end
    end
    bc.RunTimeVars.ParameterNames(1, pwvar) = pn(1, pwvar);
end

% update object
bc.Cond = ne_methods.catstruct(bc.Cond(:)', newc(:));
bc.NrOfConditions = numel(bc.Cond);
xo.C = bc;
