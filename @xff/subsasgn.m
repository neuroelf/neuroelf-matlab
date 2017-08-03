function xo = subsasgn(xo, S, V)
%SUBSASGN  Subscript-based assignment for XFF objects.
%   OBJ = SUBSASGN(OBJ, S, V) is being called for constructs such as
%
%   OBJ.Field = V; % S = struct('type', '.', 'subs', 'Field')
%
%   Further subscript referencing is allowed and will be passed on with
%   the generic xsubsasgn function (to allow objects within objects).
%
%   Calling this function manually is not advised.
%
%   See also XFF, XSUBSASGN

% Version:  v1.1
% Build:    16022911
% Date:     Feb-29 2016, 11:49 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2014, 2016, Jochen Weber
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

% global factory
global xffsngl;

% if this was called in mistake (e.g. 3rd argument is xff)
if ~isa(xo, 'xff')
    try
        xo = builtin('subsasgn', xo, S, V);
        return;
    catch xfferror
        rethrow(xfferror);
    end
end

% decide on kind of subscripting
slen  = numel(S);
stype = S(1).type;
ssubs = S(1).subs;
switch (stype)

    % typical struct
    case '.'

        % allow cell with one char call type
        if iscell(ssubs) && numel(ssubs) == 1
            ssubs = ssubs{1};
        end

        % only works for singular object
        if numel(xo) ~= 1 || xo.L(1) == 'X'
            error('neuroelf:xff:badObject', ...
                'Struct notation only works on singular non-ROOT objects.');
        end

        % everything BUT Handles (sub-set)
        if ~strcmp(ssubs, 'Handles')

            % only allow already existing fields -> direct content
            fgood = fieldnames(xo.C);
            ffound = strcmpi(ssubs, fgood);

            % if not found
            if ~any(ffound)

                % also test RunTimeVars content
                fgood = fieldnames(xo.C.RunTimeVars);
                ffound = strcmpi(ssubs, fgood);

                % bail out if still not found
                if ~any(ffound)
                    error('neuroelf:xff:invalidProperty', ...
                        'Non-existing property for this xff type.');
                end

                % get correct case fieldname
                ffound = find(ffound);
                ssubs = fgood{ffound(1)};

                % copied implementation (to avoid further if/else/end's!)
                % set complete property
                if slen == 1

                    % try setting value
                    try
                        if ~isequal(size(V), [0, 0]) || ~isa(V, 'double')
                            xo.C.RunTimeVars.(ssubs) = V;
                        else
                            xo.C.RunTimeVars.(ssubs) = [];
                        end
                        xo.H.RunTimeVarsSaved = false;
                    catch xfferror
                        rethrow(xfferror);
                    end

                % set sub value
                else

                    % try getting, altering and re-setting value
                    try
                        xo.C.RunTimeVars.(ssubs) = xsubsasgn( ...
                            xo.C.RunTimeVars.(ssubs), S(2:end), V);
                        xo.H.RunTimeVarsSaved = false;
                    catch xfferror
                        rethrow(xfferror);
                    end
                end

                % return early!
                return;

            % found in actual content
            else

                % get correct case fieldname
                ffound = find(ffound);
                ssubs = fgood{ffound(1)};
            end
        end

        % for update
        ftype = lower(xo.S.Extensions{1});
        if isfield(xffsngl.FM, ftype) && ...
            isfield(xffsngl.FM.(ftype), 'update') && ...
            isfield(xffsngl.CONF.update, ftype) && ...
            xffsngl.CONF.update.(ftype) && ...
           ~strcmp(ssubs, 'RunTimeVars') && ...
           ~strcmp(ssubs, 'Handles')
            oV = xo.C.(ssubs);
        end

        % set complete property
        if slen == 1

            % try setting value
            try
                if ~isequal(size(V), [0, 0]) || ~isa(V, 'double')
                    xo.C.(ssubs) = V;
                else
                    xo.C.(ssubs) = [];
                end
            catch xfferror
                rethrow(xfferror);
            end

        % set sub value
        else

            % try getting, altering and re-setting value
            try
                if slen ~= 2 || ~strcmp(S(2).type, '()')
                    xo.C.(ssubs) = xsubsasgn(xo.C.(ssubs), S(2:end), V);
                else
                    if ~isequal(size(V), [0, 0]) || ~isa(V, 'double')
                        xo.C.(ssubs)(S(2).subs{:}) = V;
                    else
                        xo.C.(ssubs)(S(2).subs{:}) = [];
                    end
                end
            catch xfferror

                % handle access?
                if strcmp(ssubs, 'Handles') && strcmp(S(2).type, '.')

                    % try to pass to handle
                    try
                        xo.H = xsubsasgn(xo.H, S(2:end), V);
                        return;
                    catch xfferror
                        rethrow(xfferror);
                    end

                % handle cell-array extension of string in obj.Field
                elseif isfield(xo.C, ssubs) && ...
                    strcmp(S(end).type, '{}') && ...
                    isa(V, 'char')

                    % try to use on cellstr
                    try
                        oV = subsref(xo.C, S(1:end-1));
                        if isa(oV, 'char') && size(oV, 1) == 1
                            oV = subsasgn({oV}, S(end), V);
                            if size(oV, 1) == 1 && ...
                                length(oV) == size(oV, 2)
                                oV = oV(:);
                            end
                            xo.C.(ssubs) = subsasgn( ...
                                xo.C.(ssubs), S(2:end-1), oV);
                        else
                            error('neuroelf:xff:charToCellError', ...
                                'Cannot extend single string to cell.' ...
                            );
                        end
                    catch xfferror
                        rethrow(xfferror);
                    end

                % otherwise throw error
                else
                    rethrow(xfferror);
                end
            end
        end

        % perform obj_Update call ?
        ftype = lower(xo.S.Extensions{1});
        if isfield(xffsngl.FM, ftype) && ...
            isfield(xffsngl.FM.(ftype), 'update') && ...
            isfield(xffsngl.CONF.update, ftype) && ...
            xffsngl.CONF.update.(ftype) && ...
           ~strcmp(ssubs, 'RunTimeVars')
            try
                eval([xffsngl.FM.(ftype).update{1} '(xo, ssubs, S, oV);']);
            catch xfferror
                error('neuroelf:xff:objectUpdateError', ...
                    'Error performing object update: ''%s''.', ...
                    xfferror.message);
            end
        end

        % set RunTimeVars flag
        if strcmp(ssubs, 'RunTimeVars')
            xo.H.RunTimeVarsSaved = false;
        end

    % indexing requested
    case {'()'}

        % we need non-empty cell subscript
        if ~iscell(ssubs) || isempty(ssubs)
            error('neuroelf:xff:badSubsRef', 'Can''t index into xff array.');
        end

        % try to retrieve subscripted array
        try
            subset = subsref(xo, S(1));
        catch xfferror
            error( ...
                'xff:BadSubsRef', ...
                'Invalid subscript (%s). Dynamic growing unavailable.', ...
                xfferror.message ...
            );
        end

        % no further subsasgn requested
        if slen == 1
            if isempty(V)
                xo(ssubs{:}) = [];
                return;
            end
            if ~isa(V, 'xff')
                error('neuroelf:xff:badSubsAsgnValue', ...
                    'Class mismatch error or invalid object.');
            end

            % try to assign new objects into matrix
            try
                xo(ssubs{:}) = V;
            catch xfferror
                error( ...
                    'xff:BadSubsAsgnIndex', ...
                    'Couldn''t assign partial object matrix (%s).', ...
                    xfferror.message ...
                );
            end
            return;
        end

        % further indexing only allowed for single object
        if numel(subset) ~= 1 || ~isvalid(subset) || subset.L(1) == 'X'
            error('neuroelf:xff:invalidObjectSize', ...
                'Subscripting only works for singular, non-ROOT objects.');
        end

        % try subsasgn
        try
            xo(ssubs{:}) = xsubsasgn(subset, S(2:end), V);
        catch xfferror
            error('xff:BadSubsAsgnSubs', ...
                'Error passing subsasgn to object (%s).', ...
                xfferror.message ...
            );
        end

    % generally wrong
    otherwise
        error('neuroelf:xff:badSubsAsgn', ...
            'Only struct notation allowed to set values.');
end
