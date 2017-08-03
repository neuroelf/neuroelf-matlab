function [varargout] = subsref(xo, S)
%SUBSREF  Subscript read-access method for XFF objects.
%   VALUE = SUBSREF(OBJ, S) is being called for constructs such as
%
%   VALUE = OBJ.Field; % S = struct('type', '.', 'subs', 'Field')
%
%   Further subscript referencing is allowed and will be passed on with
%   the generic xsubsref function (to allow objects within objects).
%
%   Calling this function manually is not advised.
%
%   See also XFF, XSUBSREF

% Version:  v1.1
% Build:    16032821
% Date:     Mar-28 2016, 9:03 PM EST
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

% for multiple objects
nfile = numel(xo);
if nfile > 1 && ~strcmp(S(1).type, '()')

    % prepare output (in first output!)
    varargout = cell(size(xo));

    % iterate
    for oc = 1:nfile

        % attempt
        try

            % to pass on
            varargout{oc} = xsubsref(xo(oc), S);

        % throw error
        catch xfferror
            error('neuroelf:xff:internalError', 'Error passing subsref to object %d: %s.', oc, xfferror.message);
        end
    end

    % and return
    return;
end

% default varargout
varargout = cell(1, nargout);

% number and type of first subsref
slen  = numel(S);
stype = S(1).type;
ssubs = S(1).subs;
switch (stype)

    % struct notation (object must be 1x1 here)
    case '.'

        % allow cell with one char call type
        if iscell(ssubs) && numel(ssubs) == 1
            ssubs = ssubs{1};
        end

        % get file type for methods
        ftype = lower(xo.S.Extensions{1});
        fm = xffsngl.FM;

        % try different things
        try

            % object-type specific methods (case insensitive)
            lsubs = lower(ssubs);
            if isfield(fm, ftype) && isfield(fm.(ftype), lsubs)
                tfm = fm.(ftype).(lsubs);

                % with arguments
                if slen > 1 && strcmp(S(2).type, '()')

                    % grab those
                    fargs = S(2).subs;

                    % and remove from stack (for additional subsref with .)
                    S(2)  = [];
                    slen  = numel(S);

                % without arguments
                else
                    fargs = {};
                end

                % more subsref-ing
                if slen > 1

                    % simply call with one output (rest couldn't be used)
                    [varargout{1}] = feval(tfm{6}, xo, fargs{:});

                % no more subs-refing
                else

                    % with outputs
                    if nargout > 0
                        [varargout{1:nargout}] = feval(tfm{6}, xo, fargs{:});

                    % force at least one output
                    else
                        [varargout{1}] = feval(tfm{6}, xo, fargs{:});
                    end

                    % and return early
                    return;
                end

            % methods of "any type"
            elseif isfield(fm.aft, lsubs) && ...
                (any(strcmpi(ftype, fm.aft.(lsubs){5})) || strcmpi(fm.aft.(lsubs){5}{1}, 'all'))

                % follow same logic as above
                tfm = fm.aft.(lsubs);
                if slen > 1 && strcmp(S(2).type, '()')
                    fargs = S(2).subs;
                    S(2)  = [];
                    slen  = length(S);
                else
                    fargs = {};
                end
                if slen > 1
                    [varargout{1}] = feval(tfm{6}, xo, fargs{:});
                else
                    if nargout > 1
                        [varargout{1:nargout}] = feval(tfm{6}, xo, fargs{:});
                    else
                        [varargout{1}] = feval(tfm{6}, xo, fargs{:});
                    end
                    return;
                end

            % regular field    
            elseif isfield(xo.C, ssubs)
                varargout{1} = xo.C.(ssubs);

            % field in RunTimeVars
            elseif isfield(xo.C.RunTimeVars, ssubs)
                varargout{1} = xo.C.RunTimeVars.(ssubs);

            % test special DefaultProperty from BFF/TFF code
            elseif numel(xo.S.DefaultProperty) > 1 && ...
                sum(xo.S.DefaultProperty{1} == '%') == 1 && ...
                sum(xo.S.DefaultProperty{2} == '%') == 1

                % generate expression
                subsrepl = eval(strrep(strrep(xo.S.DefaultProperty{2}, ...
                    '%', ssubs), '@', 'xo.C.'));

                % then evaluate
                varargout{1} = eval(strrep(strrep(xo.S.DefaultProperty{1}, ...
                    '%', 'subsrepl'), '@', 'xo.C.'));

            % try a "catch-abbreviation"
            else

                % begin with object-specific methods list
                meths = {};
                if isfield(fm, ftype)
                    meths = fieldnames(xffsngl.FM.(ftype));
                end

                % add the methods of all objects
                meths = [meths; fieldnames(xffsngl.FM.aft)];

                % try to find match
                methi = find(~cellfun('isempty', regexpi(meths, ['^', lower(ssubs)])));

                % not exactly one?
                if numel(methi) ~= 1
                    error('neuroelf:xff:unknownFieldOrMethod', ...
                        'Field or method %s unknown for type %s objects.', ssubs, upper(ftype));
                end

                % yet not available
                if ~isfield(xffsngl.FM.(ftype), meths{methi}) && ...
                   ~any(strcmpi(ftype, fm.aft.(meths{methi}){5})) && ...
                   ~strcmpi(fm.aft.(meths{methi}){5}{1}, 'all')
                    error('neuroelf:xff:methodUnavailable', ...
                        'Method %s not available for type %s objects.', ssubs, upper(ftype));
                end

                % re-call with proper content
                S(1).subs = meths{methi};
                if slen > 1
                    varargout{1} = xsubsref(xo, S);
                else
                    if nargout > 1
                        [varargout{1:nargout}] = xsubsref(xo, S);
                    else
                        varargout{1} = xsubsref(xo, S);
                    end
                end

                % this took care of everything (in this branch!)
                return;
            end

        % deal with errors
        catch xfferror
            rethrow(xfferror);
        end

        % more sub-indexing ?
        if slen > 1
            try

                % make sure that multiple arguments are OK
                if nargout > 0
                    [varargout{1:nargout}] = xsubsref(varargout{1}, S(2:end));
                else
                    varargout{1} = xsubsref(varargout{1}, S(2:end));
                end

            % deal with errors (sub-indexing)
            catch xfferror

                % still allow VALUE{1} for strings
                if ischar(varargout{1}) && numel(S) == 2 && strcmp(S(2).type, '{}')
                    try
                        varargout{1} = xsubsref(varargout(1), S(2));
                    catch xfferror
                        rethrow(xfferror);
                    end
                    return;
                end

                % generate an error message
                rethrow(xfferror);
            end
        end

    % regular subindexing
    case '()'

        % we need non-empty cell subscript
        if ~iscell(ssubs) || isempty(ssubs)
            error('neuroelf:xff:badSubsRef', 'Can''t index into xff array.');
        end

        % for singular object, try default property!
        if numel(xo) == 1
            if numel(xo.S.DefaultProperty) == 1
                try
                    Sa = S(1);
                    Sa.type = '.';
                    Sa.subs = xo.S.DefaultProperty{1};
                    S = [Sa; S(:)];
                    if nargout < 1
                        varargout{1} = xsubsref(xo, S);
                    else
                        [varargout{1:nargout}] = xsubsref(xo, S);
                    end
                    return;
                catch xfferror
                    rethrow(xfferror);
                end
            end
        end

        % try to retrieve subscripted matrix
        try
            subset = xo(ssubs{:});
        catch xfferror
            rethrow(xfferror);
        end

        % return sub-matrix if only one subscript
        if slen == 1
            varargout{1} = subset;
            return;
        end

        % try to pass subscripts
        try
            varargout = cell(size(subset));
            for oc = 1:numel(subset)
                varargout{oc} = xsubsref(subset(oc), S(2:end));
            end
        catch xfferror
            rethrow(xfferror);
        end

    % generally impermissible
    otherwise
        error('neuroelf:xff:badSubsRef', 'Only struct and array index subsref allowed.');
end
