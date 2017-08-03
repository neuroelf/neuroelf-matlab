function [sdm, t, ixx, ixxt] = poolnonsingletrial(sdm, preds, spred, patt)
% poolsingletrial  - pool all but one single trials of an SDM
%
% FORMAT:       [sdm, t, ixx, ixxt] = poolnonsingletrial(sdm, preds, spred [, patt])
%
% Input fields:
%
%       sdm         TxP (time-by-predictors) double design matrix
%       preds       predictor names
%       spred       single-trial predictor (either name or number)
%       patt        optional detection pattern, default: '_T\d+$'
%
% Output fields:
%
%       sdm         adapted design matrix (spred will be the first!)
%       t           if requested, transposed design matrix
%       ixx         if requested, inverse of covariance matrix
%       ixxt        if requested, ixx * t
%
% Note: the pattern MUST end in a '$'. Also, any all-zero column will be
%       removed!

% Version:  v0.9d
% Build:    14071115
% Date:     Jul-11 2014, 3:31 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, Jochen Weber
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
if nargin < 3 || ...
   ~isa(sdm, 'double') || ...
    ndims(sdm) ~= 2 || ...
    isempty(sdm) || ...
   ~iscell(preds) || ...
    numel(preds) ~= size(sdm, 2) || ...
   ((~isa(spred, 'double') || ...
     numel(spred) ~= 1 || ...
     isinf(spred) || ...
     isnan(spred) || ...
     spred < 1 || ...
     spred > numel(preds)) && ...
    (~ischar(spred) || ...
     isempty(spred)))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
preds = preds(:);
if ischar(spred)
    spred = find(strcmp(preds, spred(:)'));
    if numel(spred) ~= 1
        error( ...
            'neuroelf:BadArgument', ...
            'Predictor name not found.' ...
        );
    end
end
if nargin < 4 || ...
    ischar(patt) || ...
    isempty(patt) || ...
    isempty(regexpi(patt(:)', '\\d+')) || ...
    patt(end) ~= '$'
    patt = '_T\d+$';
else
    patt = patt(:)';
end

% re-order
npred = setdiff(1:numel(preds), spred);
sdm = sdm(:,[spred, npred]);
preds = preds([spred, npred], 1);

% match predictors
stp = (~cellfun('isempty', regexpi(preds, patt)));

% unique
ucs = unique(regexprep(preds(stp), patt, ''));
ucp = lower(regexprep(ucs, '([\[\]\{\}\:\.\?\!\+\*\^\$\\])', '\\$1'));

% iterate over each of those
for cc = 1:numel(ucs)

    % match
    stp = (~cellfun('isempty', regexpi(preds, [ucp{cc} patt])));

    % do not combine first condition
    stp(1) = false;

    % combine
    if sum(stp) > 1

        % find non-matches
        ktp = find(~stp);

        % combine
        sdm = [sdm(:, 1), sum(sdm(:, stp), 2), sdm(:, ktp(2:end))];
        preds = preds([1, 1, lsqueeze(ktp(2:end))'], 1);
        preds(2) = ucp(cc);
    end
end

% outputs
if nargout > 1
    t = sdm';
    if nargout > 2
        if any(sdm(:, 1) ~= 0) && ...
           ~any(isinf(sdm(:)) | isnan(sdm(:)))
            az = find(all(sdm == 0, 1));
            if any(az)
                sdm(:, az) = [];
                t = sdm';
            end
            ixx = inv(t * sdm);
        else
            ixx = zeros(size(t, 1), size(sdm, 2));
        end
        if nargout > 3
            ixxt = ixx * t;
        end
    end
end
