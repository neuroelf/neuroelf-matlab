function xo = aft_RemoveEntries(xo, eid, ur)
% AFT::RemoveEntries  - remove entrie/s from list with names
%
% FORMAT:       [obj = ] obj.RemoveEntries(eid [, userxi]);
%
% Input fields:
%
%       eid         either single entry (char) or list of entries (cell)
%       userxi      flag, use regexpi instead of strcmpi (default: false)
%
% Output fields:
%
%       obj         altered object
%
% Note: this method only works for objects with a list (struct) of items
%       having a suitable Name property (field)
%
% TYPES: PRT, SDM, SMP, VMP
%
% Using: multimatch.

% Version:  v1.1
% Build:    16020214
% Date:     Feb-02 2016, 2:28 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2015, 2016, Jochen Weber
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

% check arguments
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, {'prt', 'sdm', 'smp', 'vmp'}) || ...
   ((~ischar(eid) || isempty(eid) || numel(eid) ~= size(eid, 2)) && ...
    (~iscell(eid) || isempty(eid) || ~ischar(eid{1}) || isempty(eid{1})) && ~isa(eid, 'double'))
    error('neuroelf:xff:badArgument', 'Bad object or argument in call.');
end
if ischar(eid)
    eid = {eid};
elseif isa(eid, 'double')
    eid = eid(:)';
    eid(isinf(eid) | isnan(eid) | eid < 1) = [];
    eid = unique(floor(eid));
end
if nargin < 3 || ~islogical(ur) || numel(ur) ~= 1
    ur = false;
end

% depending on object type
bc = xo.C;
type = lower(xo.S.Extensions{1});
switch (type)

    % PRT
    case 'prt'

        % condition names must be treated specifically
        enames = {bc.Cond(:).ConditionName};
        try
            for ec = 1:numel(enames)
                enames{ec} = enames{ec}{1};
            end
        catch xfferror
            error('neuroelf:xff:badObject', ...
                'Invalid ConditionName for at least one condition: %s.', xfferror.message);
        end

        % match
        if iscell(eid)
            remid = ne_methods.multimatch(enames, eid, ur);
        else
            eid(eid > numel(enames)) = [];
            remid = zeros(numel(enames), 1);
            if ~isempty(eid)
                remid(eid) = 1:numel(eid);
            end
        end

        % and remove matches
        bc.Cond(remid > 0) = [];
        bc.NrOfConditions = numel(bc.Cond);

    % SDM
    case 'sdm'

        % match and remove (more fields to do)
        if iscell(eid)
            remid = ne_methods.multimatch(bc.PredictorNames(:), eid, ur);
        else
            eid(eid > numel(bc.PredictorNames)) = [];
            remid = zeros(numel(bc.PredictorNames), 1);
            if ~isempty(eid)
                remid(eid) = 1:numel(eid);
            end
        end
        bc.PredictorNames(remid > 0) = [];
        bc.NrOfPredictors = numel(bc.PredictorNames);
        bc.FirstConfoundPredictor = ...
            bc.FirstConfoundPredictor - sum(remid(1:bc.FirstConfoundPredictor-1) > 0);
        bc.PredictorColors(remid > 0, :) = [];
        bc.SDMMatrix(:, remid > 0) = [];
        bc.RTCMatrix(:, remid(1:size(bc.RTCMatrix, 2)) > 0) = [];

    % map formats
    case {'smp', 'vmp'}

        % match and remove (and patch NrOfMaps field)
        if iscell(eid)
            bc.Map(ne_methods.multimatch({bc.Map(:).Name}, eid, ur) > 0) = [];
        else
            eid(eid > numel(bc.Map)) = [];
            if ~isempty(eid)
                bc.Map(eid) = [];
            end
        end
        bc.NrOfMaps = numel(bc.Map);
end

% set content back
xo.C = bc;
