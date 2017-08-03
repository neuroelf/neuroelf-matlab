function xfigurecallback(cbstr, fig, cbo, setparent)
% xfigure::docallback  - make sure a callback has its own workspace
%
% Since this callback is only for internal use of the xfigure class
% it must not be used directly!
%
% See also @xfigure

% Version:  v1.1
% Build:    16041810
% Date:     Apr-18 2016, 10:12 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2016, Jochen Weber
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

% set references
try

    % get object (this)
    this = xfigure(cbo);
    if ~isxfigure(this, true)
        error('INVALID_HANDLE');
    end

    % set handles
    gcf  = fig;
    gcbf = fig;
    gcbo = cbo;

    % set parent?
    if nargin > 3 && ~isempty(setparent)
        parent = xfigure(get(gcbo, 'Parent'));
        thisfig = xfigure(gcf);
    end

    % evaluate callback
    eval(cbstr);

% deal with problems
catch ne_eo;
    warning('neuroelf:xfigure:badHandle', ...
        'Callback of handle failed: ''%s''.', ne_eo.message);
end
