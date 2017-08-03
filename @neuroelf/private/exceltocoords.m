function coords = exceltocoords(excelrange)
% exceltocoords  - converts MS-Excel notation into 1x4 double array
%
% FORMAT:       coords = exceltocoords(range)
%
% Input fields:
%
%       range       string representing an excel range
%
% Output fields:
%
%       coords      1x4 coords spec for range
%
% See also acsvread

% Version:  v0.9b
% Build:    11050301
% Date:     Apr-09 2011, 2:11 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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
   ~ischar(excelrange)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing Excel notation.' ...
    );
end

% split at colon
excelrange = splittocell(excelrange, ':', 1, 1);
if numel(excelrange) < 2
    excelrange{2} = excelrange{1};
end

% set defaults
ulx = 1;
uly = 1;
lrx = 256;
lry = 65535;

% get corner values
try
    [ulx, uly] = i_splitexcelcoords(excelrange{1}, 0);
    [lrx, lry] = i_splitexcelcoords(excelrange{2}, 1);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    coords = [ulx, uly, lrx, lry];
    return;
end

if lrx < ulx
    lrx = ulx;
end
if lry < uly
    lry = uly;
end
coords = [ulx, uly, lrx, lry];


%  internal function: convert 'B4' -> [2, 4]


function [cx, cy] = i_splitexcelcoords(excelcoord, minmax)

    % persistent multipliers
    persistent i_sec_m;
    if isempty(i_sec_m)
        i_sec_m = {[1, 26], [1, 10, 100, 1000, 10000]};
    end

    % default: don't use minmax setting for either column or row
    umm = false(1, 2);

    % get letter and number parts
    [cpr{1:3}] = regexpi(excelcoord, '^\s*([a-z]+)?(\d+)?\s*$');

    % if pattern didn't match at all
    if isempty(cpr{3}) || ...
       ~all(size(cpr{3}{1}) == 2)
        umm(:) = true;

    % pattern matched, get parts
    else
        px = double(lower(excelcoord(cpr{3}{1}(1, 1):cpr{3}{1}(1, 2)))) - 96;
        if isempty(px) || ...
            numel(px) > 2
            umm(1) = true;
        end
        py = double(excelcoord(cpr{3}{1}(2, 1):cpr{3}{1}(2, 2))) - 48;
        if isempty(py) || ...
            numel(py) > 5
            umm(2) = true;
        end
    end

    %  convert alpha to numeric
    if ~umm(1)
        cx = min(256, i_sec_m{1}(1:numel(px)) * px(end:-1:1)');
        if cx > 256
            cx = 256;
        end
    elseif minmax
        cx = 256;
    else
        cx = 1;
    end

    % convert numeral
    if ~umm(2)
        cy = min(65535, i_sec_m{2}(1:numel(py)) * py(end:-1:1)');
    elseif minmax
        cy = 65535;
    else
        cy = 1;
    end

% end of function i_splitexcelcoords
