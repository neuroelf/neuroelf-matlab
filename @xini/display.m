function display(S)
% xini::display  - overloaded method
%
% displays the content of an xini object
%
% FORMAT:       display(IniFileObject)
%
% Input fields:
%
%       IniFileObject   object of which content should be displayed

% Version:  v0.9d
% Build:    14090711
% Date:     Sep-07 2014, 11:56 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010 - 2014, Jochen Weber
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

% for sizes check
num_h = numel(S);
isf = isxini(S);

% if we have output arguments
if num_h == 1
    if ~isf
        error( ...
            'xini:InvalidObject', ...
            'Cannot display information on invalid object.' ...
        );
    end
    try
        disp(xini(S, 'display'));
    catch ne_eo;
        fprintf('Cannot display content: %s.\n', ne_eo.message);
    end
else
    for hc = 1:num_h
        if ~isf(hc)
            continue;
        end
        disp('----------------------------------');
        try
            disp(xini(S(hc), 'getfilename'));
            disp(xini(S(hc), 'display'));
        catch ne_eo;
            rethrow(ne_eo);
        end
    end
end
