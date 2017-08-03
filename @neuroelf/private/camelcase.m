function l = camelcase(l, opts)
% camelcase  - uppercase the first letter in each word
%
% FORMAT:       c = camelcase(label [, opts])
%
% Input fields:
%
%       label       single string or cell array of strings
%       opts        optional settings
%        .kpspc     keep spaces (default: true)
%        .mklabel   apply makelabel to the output (default: false)
%
% Output fields:
%
%       c           camelcase'd label(s)

% Version:  v0.9d
% Build:    14072210
% Date:     Jul-22 2014, 10:11 AM EST
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
if nargin < 1 || ...
    isempty(l) || ...
   (~ischar(l) && ...
    ~iscell(l)) || ...
   (iscell(l) && ...
    (any(~cellfun(@ischar, l(:))) || ...
     any(cellfun('isempty', l(:)))))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'kpspc') || ...
   ~islogical(opts.kpspc) || ...
    numel(opts.kpspc) ~= 1
    opts.kpspc = true;
end
if ~isfield(opts, 'mklabel') || ...
   ~islogical(opts.mklabel) || ...
    numel(opts.mklabel) ~= 1
    opts.mklabel = false;
end

% deblank
l = ddeblank(l);

% cell array
if iscell(l)

    % parse each (see below)
    for lc = 1:numel(l)
        ll = lower(l{lc}(:));
        lp = lsqueeze(find(ll >= 'a' & ll <= 'z'));
        nlp = [0; lsqueeze(find(ll < 'a' | ll > 'z'))];
        flp = intersect(nlp + 1, lp);
        ll(flp) = upper(ll(flp));
        if ~opts.kpspc
            ll(ll <= ' ' | ll == '_') = [];
        end
        if opts.mklabel
            l{lc} = makelabel(ll');
        else
            l{lc} = ll';
        end
    end
else

    % find letters and non-letters (incl. index 0)
    l = lower(l(:));
    lp = lsqueeze(find(l >= 'a' & l <= 'z'));
    nlp = [0; lsqueeze(find(l < 'a' | l > 'z'))];

    % find first letters
    flp = intersect(nlp + 1, lp);

    % uppercase those
    l(flp) = upper(l(flp));

    % don't keep spaces
    if ~opts.kpspc
        l(l <= ' ' | l == '_') = [];
    end

    % make label
    if opts.mklabel
        l = makelabel(l');
    else
        l = l';
    end
end
