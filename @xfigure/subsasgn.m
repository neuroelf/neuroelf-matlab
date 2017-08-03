function xo = subsasgn(xo, S, V)
% xfigure::subsasgn  - set properties on objects
%
% FORMAT:       FigureObject.PropertyName = Value

% Version:  v1.1
% Build:    16042321
% Date:     Apr-23 2016, 9:09 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% class check
if nargin > 2 && ~isa(xo, 'xfigure')
    try
        xo = builtin('subsasgn', xo, S, V);
    catch xfigerror
        rethrow(xfigerror);
    end
    return;
end

% argument check
if isempty(S)
    error('neuroelf:xfigure:badSubsAsgn', 'S struct may not be empty.');
end
slen  = numel(S);
stype = S(1).type;
ssubs = S(1).subs;

% decide on kind of subscripting, first struct-like notation
switch stype
    case '.'

        % allow cell with one char call type
        if iscell(ssubs) && numel(ssubs) == 1
            ssubs = ssubs{1};
        end

        % check for subscript type
        if ~ischar(ssubs) || isempty(ssubs)
            error('neuroelf:xfigure:badSubsAsgn', ...
                'Struct notation needs a non-empty char property.');
        end

        % only works for singular object
        if numel(xo) ~= 1
            error('neuroelf:xfigure:invalidObjectSize', ...
                'Struct notation only works on singular objects.');
        end

        % make content linear
        ssubs = ssubs(:)';

        % set complete property
        if slen == 1

            % try setting value
            try
                xoT = get(xo.H, 'Type');
                if strcmpi(xoT, 'axes')
                    if strcmpi(ssubs, 'visible')
                        set(get(xo.H, 'Children'), 'Visible', V);
                        setprop(xo, ssubs, V);
                    elseif any(strcmpi(ssubs, {'callback', 'enable'}))
                        setprop(xo, ssubs, V);
                    elseif strcmpi(ssubs, 'position') && strcmpi(xo.X.loadprops.xtype, 'xprogress') && ...
                        numel(V) == 4
                        V = V(:)';
                        set(xo.H, ssubs, V);
                        sah = get(xo.X.uicprops.xchildren(1), 'Parent');
                        Vo = get(sah, ssubs);
                        Vo(1:2) = V(1:2);
                        set(sah, ssubs, Vo);
                        xo.X.loadprops.Position(1:2) = V(1:2);
                    else
                        set(xo.H, ssubs, V);
                    end
                else
                    set(xo.H, ssubs, V);
                end
            catch xfigerror
                oeo = xfigerror;
                if ~isempty(strfind(lower(xfigerror.identifier), 'invalidproperty')) || ...
                   (isempty(xfigerror.identifier) && ...
                    ~isempty(regexpi(xfigerror.message, 'no ''\w+'' property in')))
                    try
                        set(get(xo.H, 'Children'), ssubs, V);
                        if strcmpi(ssubs, 'cdata')
                            set(xo.H, 'XLim', [0.5, 0.5 + size(V, 2)]);
                            set(xo.H, 'YLim', [0.5, 0.5 + size(V, 1)]);
                            if strcmpi(xo.X.loadprops.Type, 'ximage')
                                xo.X.loadprops.ImageData = V;
                            end
                        end
                        return;
                    catch xfigerror
                        neuroelf_lasterr(xfigerror);
                    end
                end
                rethrow(oeo);
            end

        % set sub value
        else

            % for TagStruct use alternative approach
            if strcmpi(ssubs, 'tagstruct')
                ts = tagstruct(xo);
                to = subsref(ts, S(2));
                subsasgn(to, S(3:end), V);
                return;
            end

            % try getting, altering and re-setting value
            try
                cV = get(xo.H, ssubs);
                cV = subsasgn(cV, S(2:end), V);
                set(xo.H, ssubs, cV);
            catch xfigerror
                rethrow(xfigerror);
            end
        end

    % indexing requested
    case '()'

        % we need non-empty cell subscript
        if ~iscell(ssubs) || isempty(ssubs)
            error('neuroelf:xfigure:badSubsRef', ...
                'Can''t index into xfigure matrix.');
        end

        % try to retrieve subscripted matrix
        try
            subset = xo(ssubs{:});
        catch xfigerror
            error('neureolf:xfigure:badSubsRef', ...
                'Invalid subscript (%s).', xfigerror.message);
        end

        % no further subsasgn requested
        if slen == 1
            if ~strcmpi(class(V), 'xfigure')
                error('neuroelf:xfigure:badSubsAsgnValue', ...
                    'Class mismatch error.');
            end

            % try to assign new objects into matrix
            try
                if ~isempty(V)
                    xo(ssubs{:}) = V;
                else
                    xo(ssubs{:}) = [];
                end
            catch xfigerror
                error('neuroelf:xfigure:badSubsAsgnIndex', ...
                    'Couldn''t assign partial object matrix (%s).', xfigerror.message);
            end
            return;
        end

        if numel(subset) ~= 1
            error('neuroelf:xfigure:invalidObjectSize', ...
                'Further subscripting only works for singular objects.');
        end

        try
            subsasgn(subset, S(2:end), V);
        catch xfigerror
            error('neuroelf:xfigure:badSubsAsgnSubs', ...
                'Error passing subsasgn to object (%s).', xfigerror.message);
        end

    otherwise
        error('neuroelf:xfigure:badSubsAsgn', ...
            'Only struct notation allowed to set values.');
end
