% FUNCTION ne_visspmrp: visualize SPM5/8 realignment parameter files
function varargout = ne_visspmrp(varargin)

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

% global variable
global ne_gcfg;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% only if not already in a callback
if ne_gcfg.c.incb || ...
    any(strcmp('visspmrp', ne_gcfg.c.blockcb))
    return;
end

% block further executions of this function
ne_gcfg.c.blockcb{end+1} = 'visspmrp';

% folder selector
subjdir = uigetdir(pwd, 'Please select a subject''s (functional) folder...');

% check folder
if isempty(subjdir) || ...
    isequal(subjdir, 0) || ...
    exist(subjdir, 'dir') ~= 7
    ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'visspmrp')) = [];
    return;
end

% call function
output = showspmrparams(subjdir);

% output?
if nargout > 0
    varargout{1} = output;
end

% unblock execution
ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'visspmrp')) = [];
