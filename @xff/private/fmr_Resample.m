function xo2 = fmr_Resample(xo, ifunc, tr)
% FMR::Resample  - resample the TR of the STC data of an FMR
%
% FORMAT:       newfmr = fmr.Resample(ifunc, tr)
%
% Input fields:
%
%       ifunc       interpolation function,
%                   'linear', 'cubic', 'lanczosN', 'nearest'
%       tr          new TR
%
% Using: flexinterpn_method.

% Version:  v1.1
% Build:    16020311
% Date:     Feb-03 2016, 11:19 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2016, Jochen Weber
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
flexinterpn_method = ne_methods.flexinterpn_method;

% argument check
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'ftc') || ...
   ~ischar(ifunc) || ~any(strcmpi(ifunc(:)', {'linear', 'cubic', 'nearest'})) || ...
   ~isa(tr, 'double') || numel(tr) ~= 1 || isinf(tr) || isnan(tr) || tr < 100
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
ifunc = lower(ifunc(:)');
tr = round(tr);

% copy object
xo2 = aft_CopyObject(xo);
bc = xo2.C;

% deal with empty FMRs
if isempty(bc.Slice)
    bc.TR = tr;
    xo2.C = bc;
    return;
end

% number of volumes
ss = size(bc.Slice(1).STCData);

% compute resampling factor
rfac = tr / bc.TR;

% update TR
bc.TR = tr;

% create sampling grid argument
if numel(bc.Slice) > 1
    res = [Inf, Inf, Inf; 1, 1, 1; 1, 1, rfac; ss(1), ss(2), ss(3) + (0.999 * rfac)];
else
    res = [Inf, Inf, Inf, Inf; 1, 1, 1, 1; 1, 1, rfac, 1; ss(1), ss(2), ss(3) + (0.999 * rfac), 1];
end

% for each slice (or space)
for sc = 1:numel(bc.Slice)

    % if transio
    if istransio(bc.Slice(sc).STCData)

        % resolve
        bc.Slice(sc).STCData = resolve(bc.Slice(sc).STCData);
    end

    % then resample
    if isinteger(bc.Slice(sc).STCData)
        bc.Slice(sc).STCData = uint16(flexinterpn_method(bc.Slice(sc).STCData, res, ifunc));
    else
        bc.Slice(sc).STCData = single(flexinterpn_method(bc.Slice(sc).STCData, res, ifunc));
    end
end

% set final field
bc.NrOfVolumes = size(bc.Slice(1).STCData, 3);

% but content into array
xo2.C = bc;
