%XSTRUCT  Reference-based (handle-derived) implementation of struct.
%   XS = XSTRUCT('field1', VALUE1, 'field2', VALUE2, ...) creates an
%   object variable with the same syntax as a STRUCT variable, with
%   the difference that if its contents (a field) is changed at any
%   place in the code of a second copy (or in a called function),
%   these changes are reflected everywhere, i.e. the contents is
%   shared among all references to the same XSTRUCT object.
%
%   GETFIELD and SETFIELD are not overloaded, as they work already,
%   and the same goes for ISEQUAL.
%
%   The additional overhead is substantial, particularly for structs
%   with many fields. For instance, for a xstruct with 64 fields, the
%   time to either write or read a field (in its entirety, not with
%   additional subsref components) is typically 10 to 20 times slower
%   as compared with a regular struct. So, from a run-time perspective,
%   this is a trade-off between functionality and "call-by-reference"
%   behavior that may or may not be needed.
%
%   Example:
%
%   x = xstruct('field1', 1, 'field2', 2);
%   y = x;
%   y.field2 = 3;
%   x.field2
%
%   ans =
%
%        3
%
%   To create an actual copy, you can use the copy constructor syntax:
%
%   x = xstruct('field1', 1, 'field2', 2);
%   y = xstruct(x);
%   y.field3 = 3;
%   x
%
%   x =
%
%       field1: 1
%       field2: 2
%
%   See also STRUCT, FIELDNAMES, ISSTRUCT, GETFIELD, SETFIELD,
%   XSUBSASGN, XSUBSREF.

% Version:  v1.1
% Build:    16031110
% Date:     Mar-11 2016, 10:03 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2016, Jochen Weber
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

classdef xstruct < handle

    % properties (internal set of fields)
    properties (Access = 'private', Hidden = true)

        % content of struct that is being accessed
        c = struct;
    end

    % methods (access)
    methods (Access = 'public', Hidden = true)

        % constructor
        function xs = xstruct(varargin)

            % attempt to set content
            try
                if nargin > 0
                    xs.c = struct(varargin{:});
                else
                    xs.c = struct;
                end
            catch ne_eo;
                rethrow(ne_eo);
            end
        end

        % isstruct (so that other code accepts the object)
        function tf = isstruct(xs)
            if ~isa(xs, 'xstruct')
                tf = false;
            else
                tf = true;
            end
        end
        function tf = isxstruct(xs)
            if ~isa(xs, 'xstruct')
                tf = false;
            else
                tf = true;
            end
        end

        % isempty
        function tf = isempty(xs)
            tf = isempty(xs.c);
        end

        % read-access
        function [varargout] = subsref(xs, s)

            % single () expression, return xstruct (which is YET A COPY!)
            if numel(s) == 1 && ...
                strcmp(s.type, '()')

                % all -> then keep content (reshapes original though!)
                if ~isempty(s.subs) && ...
                    all(cellfun(@ischar, s.subs)) && ...
                    all(strcmp(s.subs, ':'))
                    if numel(s.subs) == 1
                        xs.c = xs.c(:);
                    end
                    varargout{1} = xs;
                    return;
                end

                % return sub-portion (allowing any kind of indexing)
                varargout{1} = xstruct(subsref(xs.c, s));
                return;
            end

            % otherwise, pass on
            try
                [varargout{1:nargout}] = xsubsref(xs.c, s);
            catch ne_eo;
                rethrow(ne_eo);
            end
        end

        % write-access
        function xs = subsasgn(xs, s, v)

            % pass on
            try
                if numel(s) == 1 && ...
                    s.type(1) == '.'
                    xs.c.(s.subs) = v;
                else
                    xs.c = xsubsasgn(xs.c, s, v);
                end
            catch ne_eo;
                rethrow(ne_eo);
            end
        end

        % isfield, rmfield
        function tf = isfield(xs, fn)
            if nargin ~= 2
                error('neuroelf:xstruct:badNumberOfInputs', 'Invalid number of inputs.');
            end
            tf = isfield(xs.c, fn);
        end
        function xs = rmfield(xs, varargin)
            if nargin < 2
                error('neuroelf:xstruct:badNumberOfInputs', 'Invalid number of inputs.');
            end
            if nargin > 2
                xs.c = rmfield(xs.c, varargin);
            else
                xs.c = rmfield(xs.c, varargin{1});
            end
        end

        % size, ndims, length, end, and numel
        function [varargout] = size(xs, dim)
            if nargin == 1
                [varargout{1:nargout}] = size(xs.c);
            else
                if numel(dim) ~= 1 || ~isnumeric(dim) || ...
                    isinf(dim) || isnan(dim) || ...
                    dim < 1 || dim ~= fix(dim)
                    error('neuroelf:xstruct:dimensionMustBePositiveInteger', 'Dimension argument must be a positive integer scalar within indexing range.');
                end
                varargout{1} = size(xs.c, dim);
            end
        end
        function d = ndims(xs)
            d = ndims(xs.c);
        end
        function l = length(xs)
            l = length(xs.c);
        end
        function e = end(xs, k, n)
            e = builtin('end', xs.c, k, n);
        end
        function n = numel(xs, varargin)
            n = numel(xs.c);
        end

        % display
        function disp(xs)

            % invalid
            if ~isvalid(xs)
                fprintf('  handle to deleted xstruct\n\n');
                return;
            end

            % pass on
            if ~isempty(xs.c)
                disp(xs.c);
            else
                dsz = sprintf('%dx', size(xs.c));
                fn = fieldnames(xs.c);
                if isempty(fn)
                    fprintf('%s xstruct array with no fields.\n', dsz(1:end-1));
                else
                    fprintf('%s xstruct array with fields:\n\n', dsz(1:end-1));
                    fprintf('    %s\n', fn{:});
                    fprintf('\n');
                end
            end
        end

        % fieldnames, methods (for tab completion)
        function names = fieldnames(xs)
            names = fieldnames(xs.c);
        end
        function names = methods(xs)
            names = fieldnames(xs.c);
        end
        function names = properties(xs)
            names = fieldnames(xs.c);
        end

        % orderfields
        function [xs, perm] = orderfields(xs, s2)
            if nargin < 1 || nargin > 2
                error('neuroelf:xstruct:badNumberOfInputs', 'Invalid number of inputs.');
            end
            try
                if nargin == 1
                    xs.c = orderfields(xs.c);
                else
                    if nargout == 1
                        xs.c = orderfields(xs.c, s2);
                    else
                        [xs.c, perm] = orderfields(xs.c, s2);
                    end
                end
            catch ne_eo;
                rethrow(ne_eo);
            end
        end

        % content (for copy)
        function s = struct(xs, varargin)
            if nargin == 1
                s = xs.c;
            else
                s = builtin('struct', xs, varargin{:});
            end
        end

        % struct2cell
        function c = struct2cell(xs)
            c = struct2cell(xs.c);
        end

        % concatenation (necessarily creates a copy!)
        function xs = cat(dim, varargin)

            % inputs
            if nargin < 2
                error('neuroelf:xstruct:badNumberOfInputs', 'Invalid number of inputs.');
            end

            % two inputs
            if nargin == 2
                xs = varargin{1};
                return;
            end

            % classes of inputs
            icl = cellfun(@class, varargin, 'UniformOutput', false);
            if ~all(strcmp(icl, 'struct') | strcmp(icl, 'xstruct'))
                error('neuroelf:xstruct:BadInput', 'Invalid input type.');
            end

            % input sizes
            isz = cellfun(@size, varargin, 'UniformOutput', false);
            if numel(unique(cellfun('prodofsize', isz))) ~= 1
                error('neuroelf:xstruct:dimensionMismatch', 'Dimensions of matrices being concatenated are not consistent.');
            end
            isz = cat(1, isz{:});
            if dim > size(isz, 2)
                isz(:, end+1:dim) = 1;
            end
            ddiff = any(diff(isz) ~= 0, 1);
            ddiff(dim) = false;
            if any(ddiff)
                error('neuroelf:xstruct:dimensionMismatch', 'Dimensions of matrices being concatenated are not consistent.');
            end

            % output size
            sdiff = isz(1, :);
            sdiff(dim) = sum(isz(:, dim));

            % fieldnames
            ifn = cellfun(@fieldnames, varargin, 'UniformOutput', false);

            % all fields the same
            if numel(unique(cellfun('prodofsize', ifn))) == 1 && ...
                numel(unique(cat(1, ifn{:}))) == numel(ifn{1}) && ...
                isequal(repmat(ifn(1), size(ifn)), ifn)

                % use struct2cell -> cell2struct
                xs = cellfun(@struct2cell, varargin, 'UniformOutput', false);
                xs = cat(dim + 1, xs{:});
                xs = xstruct(cell2struct(xs, ifn{1}, 1));
                return;
            end

            % figure out the fieldnames
            [cfn, fi1, fi2] = unique(cat(1, ifn{:}), 'stable');

            % temporarily extend to hold number of fields
            dimac = repmat({':'}, 1, numel(sdiff));
            sdiff = [numel(cfn), sdiff];

            % generate cell
            xs = cell(sdiff);

            % fill with content
            tfc = 1;
            tfi = 1;
            for ic = 1:numel(varargin)
                nfn = numel(ifn{ic});
                dimac{dim} = tfi:tfi+isz(ic,dim)-1;
                xs(fi2(tfc:tfc+nfn-1), dimac{:}) = struct2cell(varargin{ic});
                tfc = tfc + numel(ifn{ic});
                tfi = tfi + isz(ic, dim);
            end

            % convert
            xs = xstruct(cell2struct(xs, cfn, 1));
        end

        % overload horzcat and vertcat with cat (copies in output)
        function c = horzcat(varargin)
            try
                c = cat(2, varargin{:});
            catch ne_eo;
                rethrow(ne_eo);
            end
        end
        function c = vertcat(varargin)
            try
                c = cat(1, varargin{:});
            catch ne_eo;
                rethrow(ne_eo);
            end
        end

        % overload repmat (creates copy!)
        function xnew = repmat(xs, varargin)
            xnew = xstruct(repmat(xs.c, varargin{:}));
        end

        % overload any2ascii with struct notation
        function asciirep = any2ascii(xs, varargin)
            asciirep = ['[xstruct(' any2ascii(xs.c, varargin{:}) ')]'];
        end
    end
end
