function ohelp = aft_Help(xo, smeth, varargin)
% AFT::Help  - method for any xff type
%
% FORMAT:       [helptext] = obj.Help([methodname, onerror]);
%
% Input fields:
%
%       methodname  if given, only the help for one method is returned
%       onerror     checks global flag and only returns help if requested
%
% Output fields:
%
%       helptext    1xN char array describing all (or one) method(s)
%
% TYPES: ALL
%
% Using: asciiread, gluetostring, gluetostringc, splittocell, splittocellc.

% Version:  v1.1
% Build:    16020214
% Date:     Feb-02 2016, 2:33 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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

% global config and library methods
global xffsngl ne_methods;

% only valid for single file
if numel(xo) ~= 1 || ~xffisobject(xo, true)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin > 2 && ~xffsngl.CONF.settings.Behavior.HelpOnError
    ohelp = '';
    return;
end

% locate methods for this type
ftype = xo.S.Extensions{1};
meth = xffsngl.FM;
if ~isfield(meth, lower(ftype))
    ohelp = ['No methods for filetype ' upper(ftype) '.'];
    return;
end

% get methods
ameth = meth.aft;
ometh = meth.(lower(ftype));

% only single method
if nargin > 1 && ischar(smeth) && ~isempty(smeth)
    smeth = regexprep(smeth(:)', '\s*\(.*$', '');
    if isfield(ameth, lower(smeth)) || isfield(ometh, lower(smeth))

        % get single methods help
        if isfield(ometh, lower(smeth))
            mfile = ometh.(lower(smeth));
        else
            mfile = ameth.(lower(smeth));
        end
        mhelp = gethelptxt(['@xff/private/' mfile{1} '.m']);
        ohelp = ne_methods.gluetostring(ne_methods.splittocell(mhelp, char([10, 13]), 1, 1), char(10));
        if nargin > 2
            ohelp = [char([10, 10]), ohelp];
        end
    else
        ohelp = ['No ' smeth ' method for filetype ' upper(ftype) '.'];
    end
    return;
end

% combine ameth and ometh
amthf = fieldnames(ameth);
omthf = fieldnames(ometh);
aodis = [zeros(numel(amthf), 1); ones(numel(omthf), 1)];
[tmthf, tns] = sort([amthf; omthf]);
aodis = aodis(tns);

% prepare output
ohelp = cell(numel(tns), 1);

% iterate over methods
fcall = [upper(ftype) '::'];
for mc = numel(tns):-1:1

    % for known methods:
    if aodis(mc) == 1
        ohelp{mc} = sprintf('%-74s / Syntax: %s', sprintf('%-27s - %s', ...
            [fcall, regexprep(ometh.(tmthf{mc}){1}, '^[^_]+_', '')], ...
            ometh.(tmthf{mc}){3}), ometh.(tmthf{mc}){4});
    elseif any(strcmpi(ftype, ameth.(tmthf{mc}){5}))
        ohelp{mc} = sprintf('%-74s / Syntax: %s', sprintf('%-27s - %s', ...
            [fcall, regexprep(ameth.(tmthf{mc}){1}, '^[^_]+_', '')], ...
            ameth.(tmthf{mc}){3}), strrep(ameth.(tmthf{mc}){4}, 'obj.', ...
            [ftype '.']));
    else
        ohelp(mc) = [];
    end
end

% glue to string
ohelp = ne_methods.gluetostring(ohelp, char(10));



% internal function needed to get help
function ht = gethelptxt(srcfile)
global ne_methods;
ht = ne_methods.asciiread([neuroelf_path '/' srcfile]);
ht(ht == 13) = [];
ht = ne_methods.splittocellc(ht, char(10));
ht(1) = [];
for lc = 1:numel(ht)
    if isempty(ht{lc}) || ht{lc}(1) ~= '%'
        ht(lc:end) = [];
        break;
    end
end
ht = ne_methods.gluetostringc(regexprep(ht, '^\% ?', ' '), char(10));
