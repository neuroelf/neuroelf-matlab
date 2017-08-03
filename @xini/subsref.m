function [rvalue] = subsref(hIniFile, S)
% xini::subsref  - retrieve information via struct notation
%
% FORMAT:       inistruct = IniFileObject.STRUCT
%        or     section   = IniFileObject.Section
%        or     setting   = IniFileObject.Section.Setting
%
%    For multi-dimensional object arrays, the typical indexing is valid.
%    However, only singular objects may be used for instant

% Version:  v1.0
% Build:    16011316
% Date:     Jan-13 2016, 4:40 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2016, Jochen Weber
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

% get method names for alias properties
global xinimeth;
if isempty(xinimeth)
    methods(hIniFile);
end

% argument check
if isempty(S)
    error( ...
        'xini:BadSubsRef', ...
        'S struct may not be empty.' ...
    );
end
slen  = length(S);
stype = S(1).type;
ssubs = S(1).subs;

% decide on kind of subscripting, first struct-like notation
switch stype, case {'.'}

    % allow cell with one char call type
    if iscell(ssubs) && ...
        numel(ssubs) == 1
        ssubs = ssubs{1};
    end

    % check for subscript type
    if ~ischar(ssubs) || ...
        isempty(ssubs)
        error( ...
            'xini:BadSubsRef', ...
            'Struct notation needs a non-empty char subscript.' ...
        );
    end

    % only works for singular object
    if numel(hIniFile) ~= 1
        error( ...
            'xini:InvalidFileID', ...
            'Struct notation only works for valid, singular objects.' ...
        );
    end

    % make content linear
    ssubs = ssubs(:)';

    % methods prevail
    mflds = fieldnames(xinimeth);
    mfound = find(strcmpi(ssubs, mflds));
    if ~isempty(mfound) && ...
        slen == 1 && ...
       (isempty(xinimeth.(mflds{mfound})) || ...
        xinimeth.(mflds{mfound}){1}(1) == '[')
        try
            rvalue = xini(hIniFile, lower(mflds{mfound}));
        catch ne_eo;
            rethrow(ne_eo);
        end
        return;
    end

    % allow rest only for "good" objects
    if hIniFile.L == 0
        error( ...
            'xini:InvalidObject', ...
            'Other calls/requests only allowed for non-root objects.' ...
        );
    end

    % complete or section ?
    if slen == 1

        % complete requested
        if strcmp(ssubs, 'STRUCT')
            rvalue = xini(hIniFile, 'getcomplete');

        % section requested
        else
            try
                rvalue = xini(hIniFile, 'getinisection', ssubs);
            catch ne_eo;
                rethrow(ne_eo);
            end
        end

    % setting or subsetting
    else

        % method call
        if ~isempty(mfound) && ...
           strcmp(S(2).type, '()')

            % check argument count
            oargs = xinimeth.(mflds{mfound});
            fargs = S(2).subs;
            if numel(fargs) > numel(oargs)
                error( ...
                    'xini:TooManyArguments', ...
                    'Too many arguments for call to %s.', ...
                    mflds{mfound} ...
                );
            end

            % check args
            argok = true;
            for ac = 1:numel(fargs)
                oarg = oargs{ac};
                if oarg(1) == '['
                    oarg(1) = [];
                end
                if isempty(regexpi(class(fargs{ac}), oarg))
                    argok = false;
                    break;
                end
            end
            if ~argok
                error( ...
                    'xini:BadArgument', ...
                    'Invalid arguments supplied in call to %s.', ...
                    mflds{mfound} ...
                );
            end
            if numel(fargs) < numel(oargs) && ...
                oargs{numel(fargs) + 1}(1) ~= '['
                error( ...
                    'xini:BadArgument', ...
                    'Missing %d. argument of type %s in call to %s.', ...
                    numel(fargs) + 1, ...
                    oargs{numel(fargs) + 1}(3:end - 1), ...
                    mflds{mfound} ...
                );
            end

            try
                rvalue = xini(hIniFile, lower(mflds{mfound}), fargs{:});
            catch ne_eo;
                rethrow(ne_eo);
            end
            return;
        end

        % get setting
        if S(2).type(1) == '.' && ...
            ischar(S(2).subs)
            rvalue = xini(hIniFile, 'getinisetting', S(1).subs, S(2).subs);
        else
            rvalue = xini(hIniFile, 'getinisection', S(1).subs);
            rvalue = subsref(rvalue, S(2));
        end

        % subsetting requested
        if length(S) > 2
            try
                rvalue = subsref(rvalue, S(3:end));
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                error( ...
                    'xini:IllegalSubsRef', ...
                    'Couldn''t pass subsref to %s.%s.', ...
                    S(1).subs, ...
                    S(2).subs ...
                );
            end
        end
    end

% indexing requested
case {'()'}

    % we need non-empty cell subscript
    if ~iscell(ssubs) || ...
        isempty(ssubs)
        error( ...
            'xini:BadSubsRef', ...
            'Can''t index into xini matrix.' ...
        );
    end

    % convert hIniFile to struct
    sIniFile = struct(hIniFile);

    % try to retrieve subscripted matrix
    try
        subset   = subsref(sIniFile, S(1));
        hIniFile = class(subset, 'xini');
    catch ne_eo;
        error( ...
            'xini:BadSubsRef', ...
            'Invalid subscript (%s).', ...
            ne_eo.message ...
        );
    end

    % return sub-matrix if only one subscript
    if slen == 1
        rvalue = hIniFile;
        return;
    end

    if numel(hIniFile) ~= 1 || ...
       ~hIniFile.L
        error( ...
            'xini:InvalidFileID', ...
            'Further subscripting only works for valid, singular objects.' ...
        );
    end

    rvalue = subsref(hIniFile, S(2:end));

otherwise
    error( ...
        'xini:BadSubsRef', ...
        'Only struct notation allowed to retrieve values.' ...
    );
end
