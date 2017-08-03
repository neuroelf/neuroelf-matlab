function [keyPressed, rt] = ptb_waitforkey(keycodes, timeout, timein)
% ptb_waitforkey  - wait for keyboard input
%
% FORMAT:       [key, rt] = ptb_waitforkey(keycodes [, timeout [, timein]])
%
% Input fields:
%
%       keycodes    1xC key codes that are valid
%       timeout     time-out in seconds
%       timein      start time (return value of GetSecs) from prior call
%
% Output fields:
%
%       key         key name of pressed key
%       rt          reaction time (relative to timein/call)
%
% Note: if timein is not given, a call to GetSecs is done at the beginning.

% Version:  v0.9c
% Build:    12051815
% Date:     May-17 2012, 1:18 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012, Jochen Weber
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
realtimein = GetSecs;
if nargin < 3 || ...
   ~isa(timein, 'double') || ...
    numel(timein) ~= 1 || ...
    isinf(timein) || ...
    isnan(timein) || ...
    timein < 0
    timein = realtimein;
end
if nargin < 1
    error( ...
        'neuroelf:MissingArgument', ...
        'Missing argument.' ...
    );
end
if ischar(keycodes)
    keycodes = {keycodes};
end
if iscell(keycodes)
    try
        for kc = 1:numel(keycodes)
            keycodes{kc} = KbName(keycodes{kc}(:)');
        end
        keycodes = cat(2, keycodes{cellfun('isclass', keycodes, 'double')});
    catch ne_eo;
        rethrow(ne_eo);
    end
end
if ~isa(keycodes, 'double') || ...
    isempty(keycodes) || ...
    any(isinf(keycodes(:)) | isnan(keycodes(:)) | ...
        keycodes(:) < 1 | keycodes(:) > 255 | keycodes(:) ~= fix(keycodes(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2 || ...
   ~isa(timeout, 'double') || ...
    numel(timeout) ~= 1 || ...
    isnan(timeout) || ...
    timeout <= 0
    timeout = Inf;
end
timeout = timeout + timein;
keycodes = unique(keycodes(:));

% preparation
olddis = [];
[keyIsDown, secs, keyCode] = KbCheck;
if keyIsDown
    keyCode(keycodes) = 0;
    if any(keyCode)
        olddis = DisableKeysForKbCheck(find(keyCode));
        keyCode(olddis) = 1;
        DisableKeysForKbCheck(find(keyCode));
    end
end

% default
keyPressed = '';
rt = NaN;

% loop
while secs <= timeout

    % inspect keyboard
    [keyIsDown, secs, nkeyCode, deltaSecs] = KbCheck;

    % key is pressed
    if keyIsDown

        % see if a key of requested code is pressed
        keyPressedTest = keycodes(nkeyCode(keycodes) > 0);
        if ~isempty(keyPressedTest)

            % then get reaction time and out of the loop
            keyPressed = keyPressedTest;
            rt = secs - (realtimein + 0.5 * deltaSecs);
            break;
        end
    end
end

% restore disabled keys
DisableKeysForKbCheck(olddis);

% translate
if ~isempty(keyPressed)
    keyPressed = KbName(keyPressed);
end
