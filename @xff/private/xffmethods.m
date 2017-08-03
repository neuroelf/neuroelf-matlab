function ff_methods = xffmethods(ext)
% xff::_methods  - searching for methods for files
%
% THIS FUNCTION IS AN INTERNAL FUNCTION OF THE CLASS
%
% @xff
%
% AND IT SHOULD NEVER BE CALLED FROM OUTSIDE THE CLASS CODE

% Version:  v1.1
% Build:    16021216
% Date:     Feb-12 2016, 4:32 PM EST
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

% global access to methods
global ne_methods;

% create struct
ff_methods = struct;

% get path of this m-file
spfile = mfilename('fullpath');
sppath = fileparts(spfile);

% get extension list and add special type for AllFileTypes
extlist = fieldnames(ext);
extlist{end + 1} = 'aft';

% for each extension, look for methods
for ec = 1:numel(extlist)

    % get extension and build empty sub struct
    typext  = lower(extlist{ec});
    typmst  = struct;

    % lookup methods
    typmeth = dir([sppath filesep typext '_*.m']);

    % iterate over found entries
    for mc = 1:length(typmeth)

        % read file
        fcont = ne_methods.asciiread([sppath filesep typmeth(mc).name]);
        if ~isempty(strfind(fcont, char([13, 10])))
            fcont = ne_methods.splittocellc(fcont, char([13, 10]), false, false);
        else
            fcont = ne_methods.splittocellc(fcont, char(10), false, false);
        end

        % get first comment line text
        try
            fcomm = regexprep(regexprep(fcont{2}, '^.*\s+\-\s+', ''), '\s+$', '');
        catch xfferror
            neuroelf_lasterr(xfferror);
            fcomm = '';
        end

        % get FORMAT line
        fline = ne_methods.grep(fcont, ['FOR' 'MAT']);
        if ~iscell(fline)
            fline = cellstr(fline);
        end

        % get TYPES line
        if ec < numel(extlist)
            tline = {upper(typext)};
        else
            tline = ne_methods.grep(fcont, ['TYP' 'ES:']);
            if ~isempty(tline)
                tline = regexprep(tline{1}, ['^.*TYP' 'ES\:\s*'], '');
                tline(tline == ' ') = [];
                tline = ne_methods.splittocellc(tline, ',');
            end
        end

        % put method name (without extension and typext into substruct
        mname = typmeth(mc).name(1:end-2);
        uscpos = find(mname == '_');
        if ~isempty(fline)
            fargs = regexprep(fline{1}, '^[^\(]+(\([^\)]*\))?[^\)]*$', '$1');
            fline = regexprep(fline{1}, ['^.*FOR' 'MAT\:?\s*'], '');
            typmst.(lower(mname(uscpos(1)+1:end))) = ...
                {mname, fargs, fcomm, fline, tline, eval(['@' mname])};
        else
            typmst.(lower(mname(uscpos(1)+1:end))) = {mname, '(...)', fcomm, ...
                [typext '.' regexprep(typmeth(mc).name, '^\w+_', '') '(...);'], ...
                tline, eval(['@' mname])};
        end
    end

    % put into global struct
    if numel(fieldnames(typmst)) > 0
        ff_methods.(typext) = typmst;
    end
end
