function hIniFile = subsasgn(hIniFile, S, V)
% xini::subsasgn  - write settings to an ini-file
%
% FORMAT:       IniFileObject.STRUCT = inistruct;
%        or     IniFileObject.Section  = section;
%        or     IniFileObject.Section.Setting = value;
%
%    For multi-dimensional object arrays, the typical indexing is valid.
%    However, only singular objects may be used for write access.

% Version:  v1.1
% Build:    16021815
% Date:     Feb-18 2016, 3:46 PM EST
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

% class check
if nargin > 2 && ...
    ~isa(hIniFile, 'xini')
    try
        hIniFile = builtin('subsasgn', hIniFile, S, V);
    catch ne_eo;
        rethrow(ne_eo);
    end
    return;
end

% argument check
if isempty(S)
    error( ...
        'xini:BadSubsAsgn', ...
        'S struct may not be empty.' ...
    );
end
slen= length(S);
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
            'xini:BadSubsAsgn', ...
            'Struct notation needs a non-empty char subscript.' ...
        );
    end

    % only works for singular object
    if numel(hIniFile) ~= 1 || ...
       ~hIniFile.L
        error( ...
            'xini:InvalidFileID', ...
            'Struct notation only works for valid, singular objects.' ...
        );
    end

    % make content linear
    ssubs = ssubs(:)';

    % complete or section ?
    if slen == 1

        % complete requested
        if strcmp(ssubs, 'STRUCT')
            try
                xini(hIniFile, 'setinicomplete', V);
            catch ne_eo;
                rethrow(ne_eo);
            end

        % section requested
        else
            try
                xini(hIniFile, 'setinisection', S(1).subs, V);
            catch ne_eo;
                rethrow(ne_eo);
            end
        end

    % setting or subsetting
    else

        % subsetting requested
        if length(S) > 2
            try

                % 1x1 struct (without sub-indexing)
                if S(2).type(1) == '.' && ...
                    ischar(S(2).subs)

                    % get setting
                    rvalue = xini( ...
                        hIniFile, 'getinisetting', S(1).subs, S(2).subs);

                    % assign changes
                    rvalue = subsasgn(rvalue, S(3:end), V);

                    % store changes setting
                    xini( ...
                        hIniFile, 'setinisetting', S(1).subs, S(2).subs, rvalue);

                % sub-indexing
                else
                    rvalue = xini(hIniFile, 'getinisection', S(1).subs);
                    if ~isequal(size(V), [0, 0]) || ~isa(V, 'double')
                        rvalue = subsasgn(rvalue, S(2:end), V);
                    else
                        rvalue = subsasgn(rvalue, S(2:end), []);
                    end
                    xini(hIniFile, 'setinisection', S(1).subs, rvalue);
                end
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                error( ...
                    'xini:IllegalSubsAsgn', ...
                    'Couldn''t pass subsasgn to %s.%s.', ...
                    S(1).subs, ...
                    S(2).subs ...
                );
            end

        % set value directly
        else
            try
                if S(2).type(1) == '.' && ...
                    ischar(S(2).subs)
                    xini( ...
                        hIniFile, 'setinisetting', S(1).subs, S(2).subs, V);
                else
                    rvalue = xini(hIniFile, 'getinisection', S(1).subs);
                    if ~isequal(size(V), [0, 0]) || ~isa(V, 'double')
                        rvalue = subsasgn(rvalue, S(2), V);
                    else
                        rvalue = subsasgn(rvalue, S(2), []);
                    end
                    xini(hIniFile, 'setinisection', S(1).subs, rvalue);
                end
            catch ne_eo;
                rethrow(ne_eo);
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
        subset = subsref(sIniFile, S(1));
    catch ne_eo;
        error( ...
            'xini:BadSubsRef', ...
            'Invalid subscript (%s).', ...
            ne_eo.message ...
        );
    end

    % no further subsasgn requested
    if slen == 1
        if ~isa(V, 'xini')
            error( ...
                'xini:BadSubsAsgnValue', ...
                'Class mismatch error.' ...
            );
        end

        % try to assign new objects into matrix
        try
            sIniFile = subsasgn(sIniFile, S(1), struct(V));
            hIniFile = class(sIniFile, 'xini');
        catch ne_eo;
            error( ...
                'xini:BadSubsAsgnIndex', ...
                'Couldn''t assign partial object matrix (%s).', ...
                ne_eo.message ...
            );
        end
        return;
    end

    if numel(subset) ~= 1 || ...
       ~subset.L
        error( ...
            'xini:InvalidFileID', ...
            'Further subscripting only works for valid, singular objects.' ...
        );
    end

    try
        subsasgn(class(subset, 'xini'), S(2:end), V);
    catch ne_eo;
        rethrow(ne_eo);
    end

otherwise
    error( ...
        'xini:BadSubsAsgn', ...
        'Only struct notation allowed to set values.' ...
    );
end
