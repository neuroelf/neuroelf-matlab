% FUNCTION ne_echo: show function call in prompt
function ne_echo(o, c, varargin)

% Version:  v0.9b
% Build:    10081109
% Date:     Aug-11 2010, 9:00 AM EST
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

% general call
if nargin < 2

    % still sprintf
    if iscell(o)
        disp(sprintf(o{1}, o{2:end}));

    % just plain text
    else
        disp(o);
    end
% object call without arguments
elseif nargin < 3

    % just display the call
    disp(sprintf('%s.%s;', o, c));

% otherwise
else

    % print with arguments
    if nargin == 3
        disp(sprintf('%s.%s(%s);', o, c, any2ascii(varargin{1})));
    else
        a2a = cell(1, nargin - 2);
        for a2c = 1:(nargin - 2)
            if ~isxff(varargin{a2c}) || ...
               ~isempty(varargin{a2c}.FilenameOnDisk)
                a2a{a2c} = any2ascii(varargin{a2c});
            else
                a2a{a2c} = ...
                    sprintf('xff(''<%s>'')', upper(varargin{a2c}.Filetype));
            end
        end
        disp(sprintf('%s.%s(%s);', o, c, gluetostringc(a2a, ', ')));
    end
end
