% FUNCTION ne_cm_autocons: automatically add contrasts
function ne_cm_autocons(varargin)

% Version:  v1.0
% Build:    16040515
% Date:     Apr-05 2016, 3:50 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2016, Jochen Weber
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
ch = ne_gcfg.h.CM.h;
prednames = ch.Predictors.String;

% copy (used to set predictor names back in the end)
predorig = prednames;

% replace - and _ with space
prednames = strrep(strrep(prednames, '-', ' '), '_', ' ');

% make copy (used to temporarily set in predictor names)
predcopy = prednames;

% replace _T\d+ with _T
prednames = unique(regexprep(prednames, '_T\d+$', '_T'));

% now replace 'SomeText' with 'some text'
prednames = ddeblank(lower(regexprep(prednames, '(\s*|[a-z0-9_])([A-Z])', '$1 $2')));

% split into components
predparts = cellfun(@splittocell, prednames, repmat({' '}, size(prednames)), 'UniformOutput', false);

% get the number of parts
numparts = cellfun('prodofsize', predparts);

% get the median number (assuming these are the conditions we care about)
mednum = median(numparts);

% get all the predictors that match this number
matched = (numparts == mednum);
prednames = prednames(matched);
predparts = cat(1, predparts{matched});

% any parts the same?
for pc = mednum:-1:1
    if all(strcmp(predparts{1, pc}, predparts(:, pc)))
        delpart = predparts{1, pc};

        % delete from predictor names if not final t (single-trials)
        if pc < mednum || ~strcmpi(delpart, 't')
            predparts(:, pc) = [];
            prednames = regexprep(prednames, ['(\s*)' lower(delpart) '\s*'], '$1');
            mednum = mednum - 1;
        end
    end
end

% camelcase predictor names
prednames = camelcase(prednames);

% part options
partopts = cell(1, mednum);
partall = cell(1, mednum);
for pc = 1:mednum
    partopts{pc} = unique(predparts(:, pc));
    partall{pc} = sprintf('%s|', partopts{pc}{:});
    partall{pc} = ['.*(' partall{pc}(1:end-1) ').*'];
    for ppc = 1:numel(partopts{pc})
        partopts{pc}{ppc} = ['.*' partopts{pc}{ppc} '.*'];
    end
end

% total number of contrasts is the product of the number of parts
totalcons = cellfun('prodofsize', partopts);
totalcons = 1 + sum(0.5 .* totalcons .* (totalcons - 1));
connames = cell(totalcons, 1);
conform = cell(totalcons, 1);

% begin with main effect
connames{1} = 'Main Effect (all > baseline)';
conform{1} = ['^' strrep(sprintf('%s', partall{:}), '.*.*', '.*') '$'];

% remove emptys (just in case)
emptycons = cellfun('isempty', connames);
connames(emptycons) = [];
conform(emptycons) = [];

% ask
allcons = sprintf('%s\n', connames{:});
vans = questdlg(sprintf('No contrasts in current GLM. Create the following?\n\n%s', allcons), ...
    'NeuroElf - user input', 'Yes', 'No', 'Yes');
