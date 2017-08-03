function [varargout] = subsref(xo, S)
% xfigure::subsref  - retrieve property via struct notation
%
% FORMAT:       propvalue = FigureObject.PropertyName
%
% Also, the subsref construct can be used as an alternative way of
% calling an object method:
%
% FORMAT:       FigureObject.MethodName([Arguments]);
%
% Since method names are checked first, Parent() returns the parent
% object reference, not MABLAB's parent GUI handle !

% Version:  v1.1
% Build:    16042313
% Date:     Apr-23 2016, 1:25 PM EST
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
global xfiguremeth;
if isempty(xfiguremeth)
    try
        methods(xo);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
end

% output
if nargout
    varargout = cell(1, nargout);
end

% argument check
if isempty(S)
    error('neuroelf:xfigure:badSubsRef', 'S struct may not be empty.');
end
slen  = length(S);
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
            error('neuroelf:xfigure:badSubsRef', ...
                'Struct notation needs a non-empty char property.');
        end

        % only works for singular object
        if numel(xo) ~= 1
            error('xfigure:InvalidInputSize', ...
                'Struct notation works only on singular objects.');
        end

        % make content linear
        ssubs = ssubs(:)';

        % try to retrieve value
        try
            if any(strcmpi(ssubs, fieldnames(xfiguremeth.m)))
                if slen > 1 && strcmp(S(2).type, '()')
                    fargs = S(2).subs;
                    S(2)  = [];
                    slen  = length(S);
                else
                    fargs = {};
                end
                xmeth = xfiguremeth.m.(lower(ssubs));
                if xmeth{2} > 0 && nargout
                    [varargout{:}] = feval(xmeth{1}, xo, fargs{:});
                else
                    if ~strcmpi(ssubs, 'delete')
                        feval(xmeth{1}, xo, fargs{:});
                    else
                        deletefigure(xo);
                    end
                    if nargout > 0
                        varargout{1} = [];
                    end
                end
            else
                varargout{1} = getprop(xo, ssubs);
            end
        catch xfigerror
            rethrow(xfigerror);
        end

        % more sub-indexing ?
        if slen > 1
            try
                if slen > 2 && ...
                    strcmpi(ssubs, 'tagstruct') && ...
                    strcmp(S(2).type, '.') && ...
                    ischar(S(2).subs)
                    varargout{1} = varargout{1}.(S(2).subs);
                    varargout{1} = subsref(varargout{1}, S(3:end));
                else
                    varargout{1} = subsref(varargout{1}, S(2:end));
                end
            catch ne_eo;
                error( ...
                    'xfigure:IllegalSubsRef', ...
                    'Couldn''t pass further subscripting to property value: %s.', ...
                    ne_eo.message ...
                );
            end
        end

    % indexing requested
    case {'()'}

        % we need non-empty cell subscript
        if ~iscell(ssubs) || ...
            isempty(ssubs)
            error( ...
                'xfigure:BadSubsRef', ...
                'Can''t index into xfigure matrix.' ...
            );
        end

        % try to retrieve subscripted matrix
        try
            subset = xo(ssubs{:});
        catch xfigerror
            error('neuroelf:xfigure:badSubsRef', ...
                'Invalid subscript error (%s).', xfigerror.message);
        end

        % return sub-matrix if only one subscript
        if slen == 1
            varargout{1} = subset;
            return;
        end

        if numel(xo) ~= 1
            error( ...
                'xfigure:InvalidObjectSize', ...
                'Further subscripting is only valid for singular objects.' ...
            );
        end

        % try to pass subscripts
        try
            varargout{1} = subsref(subset, S(2:end));
        catch ne_eo;
            rethrow(ne_eo);
        end

    otherwise
        error( ...
            'xfigure:BadSubsRef', ...
            'Only struct notation allowed to retrieve values.' ...
        );
end
