function memusage
% memusage  - give the current memory usage of the calling function
%
% FORMAT:       memusage
%
% No input/output fields.

% Version:  v0.9d
% Build:    14090214
% Date:     Sep-02 2014, 2:04 PM EST
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

% requires no argument check

% get dbstack to get the function/line of caller, allowing for @neuroelf
d = dbstack(1);
if ~isempty(d) && ...
    strcmp(d(1).file, 'subsref.m')
    disp('memusage must be used as a function handle, not with @neuroelf.');
    return;
end

% get workspace from caller
w = evalin('caller', 'whos');

% if call not from base (prompt/UI)
if ~isempty(d)
    fprintf('%s (%d): %d bytes\n', ...
        d(1).file, d(1).line, sum(cat(1, w(:).bytes)));

% for calls from base (prompt/UI)
else
    fprintf('BASE: %d bytes\n', sum(cat(1, w(:).bytes)));
end

% make sure the screen is updated to draw text line
drawnow;
