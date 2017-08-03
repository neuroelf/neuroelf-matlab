function pv = makelabel(pv, withdots)
% makelabel  - returns a valid label from input
%
% when a valid label/identifier name is required from a string,
% this function returns a validated string that can be used.
%
% FORMAT:       vlabel = makelabel(testlabel [,withdots])
%
% Input fields:
%
%       testlabel   string to use for label generation
%       withdots    if argument is specified, struct-like
%                   identifier are accepted
%
% Output fields:
%
%       vlabel      string that can be used for a label name
%
% See also isrealvarname.

% Version:  v0.9d
% Build:    14052916
% Date:     May-29 2014, 4:28 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, Jochen Weber
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
if nargin < 1
    error( ...
        'neuroelf:TooFewArguments',...
        'Too few arguments. Try ''help %s''.',...
        mfilename ...
    );
end

% for cell arrays
if iscell(pv)

    % depending on input arguments
    if nargin == 1

        % run for each in turn
        for ec = 1:numel(pv)
            pv{ec} = makelabel(pv{ec});
        end
    else
        for ec = 1:numel(pv)
            pv{ec} = makelabel(pv{ec}, true);
        end
    end

    % return early
    return;
end

% some first checks
if ~ischar(pv)
    pv = '';
end
pv = deblank(pv(:)');

% empty variable ...
if isempty(pv)
    pv = sprintf('V_%010.0f', floor(2 ^ 31 * rand(1)));

% not with dots
elseif nargin == 1

    % return name if already is varname
    if isrealvarname(pv)
        return;
    end

    % otherwise check first character and prepend 'V_' if necessary
    pa = pv(1);
    if pa < 65 || ...
       (pa > 90 && pa < 97) || ...
        pa > 122
        pv = ['V_' pv];
    end

    % replace all bad chars with underscores
    pv(pv < 48 | (pv > 57 & pv < 65) | (pv > 90 & pv < 95) | pv == 96 | pv > 122) = '_';

    % replace double underscores with single ones
    while ~isempty(strfind(pv, '__'))
        pv = strrep(pv, '__', '_');
    end

    % get correct variable size
    pv = pv(1:min(namelengthmax, size(pv,2)));

% with dots
elseif ~isempty(withdots)
    % split to single names (at dots)
    [pcs, pcc] = splittocell(pv, '.', 1);

    % check vor every part
    for pcc = 1:pcc
        pcp = pcs{pcc};

        % allow sub-indexing (remove from name)
        pcb = regexp(pcp,'\(\d+[\,\d]*\)$');
        if ~isempty(pcb)
            pcx = pcp(pcb(1):end);
            pcp(pcb(1):end) = [];
        else
            pcx = '';
        end

        % test remaining string and recode if needed
        if ~isrealvarname(pcp)
            pcs{pcc} = [makelabel(pcp) pcx];
        end
    end

    % make a dotted expression again
    pv = gluetostring(pcs, '.');
end
