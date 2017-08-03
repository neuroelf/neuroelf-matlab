function uzstring = unzerostring(uzstring, lastzero)
% unzerostring  - removed string part after first (or last) zero
%
% FORMAT:       uzstring = unzerostring(uzstring, lastzero)
%
% Input fields:
%
%       uzstring    string to be zero-char free
%       lastzero    used last zeros instead
%
% Output fields:
%
%       uzstring    zero-char removed string

% Version:  v0.9a
% Build:    10051716
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
   ~ischar(uzstring)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid argument given to %s.', ...
        mfilename ...
    );
end
if isempty(uzstring)
    uzstring = '';
    return;
end
if nargin < 2 || ...
    isempty(lastzero) || ...
   (~isnumeric(lastzero) && ~islogical(lastzero)) || ...
   ~lastzero(1)
    lastzero = false;
else
    lastzero = true;
end
uzstring  = uzstring(:)';
duzstring = double(uzstring);

% for last zero
if lastzero

    % get length and find first non-zero char from end
    lzp = findfirst(duzstring ~= 0, -1);
    if isempty(lzp)
        uzstring = '';
    elseif lzp < numel(duzstring)
        uzstring(lzp+1:end) = [];
    end

% for first zero
else

    % find zeros
    fzp = findfirst(duzstring == 0);

    % return on no zeros
    if isempty(fzp)
        return;

    % return empty string if first char is zero
    elseif fzp == 1
        uzstring = '';
        return;
    end

    % otherwise cut string
    uzstring = uzstring(1:fzp-1);
end
