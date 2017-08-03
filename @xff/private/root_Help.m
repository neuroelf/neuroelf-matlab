function ohelp = root_Help(xo, m)
% ROOT::Help  - get Help on methods
%
% FORMAT:       [helptext] = xff.Help([typmeth]);
%
% Input fields:
%
%       typmeth     optional type or method
%
% Output fields:
%
%       help        complete help over all methods
%
% Using: gluetostringc.

% Version:  v1.1
% Build:    16021216
% Date:     Feb-12 2016, 4:27 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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
global ne_methods xffsngl;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'root')
    error('neuroelf:xff:BadArgument', 'Invalid call to ''%s''.', mfilename);
end

% get list of extensions
x = xffsngl.EXT;
xf = fieldnames(x);

% with specifics
if nargin > 1 && ischar(m) && ~isempty(m)
    m = lower(m(:)');

    % rethrow error if necessary
    try
        obj = cell(1, 1);
        
        % test for extension
        if strcmpi(m, 'root')
            ohelp = aft_Help(xffsngl.OBJS{1, 4});
            
        elseif isfield(x, m)

            % pass good object
            obj{1} = xff(['new:' m]);
            ohelp = aft_Help(obj{1});
            clearxffobjects(obj);

        % or pass to aft_Help for root
        else
            ohelp = aft_Help(xo, m);
        end
    catch xfferror;
        clearxffobjects(obj);
        rethrow(xfferror);
    end
    return;
end

% remove double entries
for xc = numel(xf):-1:2
    if x.(xf{xc}){2} == x.(xf{xc-1}){2}
        xf(xc) = [];
    end
end
xf = sort(xf);

% create new objects and get help
ohelp = cell(numel(xf), 1);
for xc = numel(xf):-1:1
    obj = [];
    try
        if strcmpi(xf{xc}, 'root')
            ohelp{xc} = aft_Help(xo);
        else
            obj = xff(['new:' xf{xc}]);
            ohelp{xc} = aft_Help(obj);
        end
    catch xfferror
        neuroelf_lasterr(xfferror);
        ohelp{xc} = '';
    end
    if ~isempty(obj)
        delete(obj);
    end
    ohelp{xc} = regexprep(ohelp{xc}, '^No methods.*$', '');
    if isempty(ohelp{xc})
        ohelp(xc) = [];
    end
end
ohelp = ne_methods.gluetostringc(ohelp, char(10));
