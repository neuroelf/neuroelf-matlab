function m = dispslicemovie(voldata, sdir, pauseval)
% dispslicemovie  - display slice movie
%
% FORMAT:       dispslicemovie(voldata, [dir [, pauseval]])
%
% Input fields:
%
%       voldata     3-D dataset
%       dir         either of [-3, -2, -1, 1, 2, 3]
%       pauseval    value to pause between two slices (default: 0.04)
%
% See also dispslice

% Version:  v0.9a
% Build:    10100115
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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
    isempty(voldata)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing input argument provided.' ...
    );
end

% resolve transio if necessary
if istransio(voldata)
    voldata = resolve(voldata);
end

% get size and necessary values
vsz = size(voldata);
if numel(vsz) < 3
    vsz(3) = 1;
end
mmm = minmaxmean(voldata, 4);
vmn = mmm(1);
vmx = mmm(2);

% no/bad direction/pauseval given
if nargin < 2 || ...
   ~isa(sdir, 'double') || ...
    numel(sdir) ~= 1 || ...
   ~any(sdir == [-4:-1, 1:4]) || ...
    abs(sdir) > numel(vsz)
    sdir = 3;
end
if nargin < 3 || ...
   ~isa(pauseval, 'double') || ...
    numel(pauseval) ~= 1 || ...
    isinf(pauseval) || ...
    isnan(pauseval) || ...
    pauseval < 0 || ...
    pauseval > 2
    pauseval = 0.04;
end

% for 4-D data, reshape to a good set
if abs(sdir) == 4
    voldata = reshape(permute(voldata, [1, 2, 4, 3]), [vsz(1:2), vsz(3) * vsz(4)]);
    sdir = 3 * sign(sdir);
elseif numel(vsz) > 3
    voldata = reshape(voldata, [vsz(1:2), vsz(3) * vsz(4)]);
end
vsz = size(voldata);

% delete dispslice figure
dispslice closefig;

% slice preselection is empty
slc = [];

% more than single direction given
if numel(sdir) > 1
    slc = sdir(2:end);
    sdir = sdir(1);
end

% reject invalid dir's
if ~any([-3:-1, 1:3] == sdir)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid dir argument.' ...
    );
end

% if negative go backwards
reversed = false;
if sdir < 0
    sdir = -sdir;
    reversed = true;
end

% force good slices
slc = intersect(slc, 1:vsz(sdir));

% preselection empty ?
if isempty(slc)

    % all slices
    slc = 1:vsz(sdir);
end

% reversed
if reversed
    slc = slc(end:-1:1);
end

% get slice size
switch (sdir)
    case {1}
        ssz = vsz(2:3);
    case {2}
        ssz = vsz([1, 3]);
    otherwise
        ssz = vsz(1:2);
end
ssx = vsz(sdir);

% create empty dispslice
dispslice(zeros(ssz));

% get image handle
ih = dispslice('image');

% required axes handle
if nargout > 0
    ia = get(ih, 'Parent');
end

% display requested slices
for sl = slc
    switch (sdir)
        case {1}
            dispslice(voldata(sl, :, :), ih, [vmn, vmx]);
        case {2}
            dispslice(voldata(:, sl, :), ih, [vmn, vmx]);
        otherwise
            dispslice(voldata(:, :, sl), ih, [vmn, vmx]);
    end
    dispslice(sprintf('Dispslicemovie: slice %d/%d', sl, ssx));
    pause(pauseval);

    % create movie ?
    if nargout > 0
        if sl == slc(1)
            m = getframe(ia);
            m(2:numel(slc)) = m;
            mi = 2;
        else
            m(mi) = getframe(ia);
            mi = mi + 1;
        end
    end
end
