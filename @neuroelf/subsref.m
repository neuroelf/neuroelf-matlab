function [varargout] = subsref(n, S)
% neuroelf::subsref  - make method (function) call to NeuroElf function
%
% FORMAT:       [varargout = ] NEUROELF_OBJECT.METHOD(varargin)
%
% METHODS:      to see a list of available methods, please use the syntax
%
%       ne_object = neuroelf % without semicolon
%
% Input/output arguments: please see help as in
%
%       ne_object.help('methodname')

% Version:  v1.1
% Build:    16041111
% Date:     Apr-11 2016, 11:04 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2014, 2016, Jochen Weber
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

% valid?
if nargin ~= 2 || ~isstruct(S) || isempty(S) || numel(S) > 2 || ...
   ~isfield(S, 'type') || ~isfield(S, 'subs') || ~strcmp(S(1).type, '.') || ...
   ~ischar(S(1).subs) || isempty(S(1).subs) || (numel(S) > 1 && ~strcmp(S(2).type, '()'))
    error('neuroelf:general:badArgument', ...
        'Bad or missing argument to neuroelf::subsref method call.');
end

% prepare outputs
varargout = cell(1, nargout);

% get struct of object
neos = struct(n);

% help?
if strcmpi(S(1).subs, 'help')

    % main help
    if numel(S) == 1
        h = asciiread([neuroelf_path filesep '@neuroelf' filesep 'neuroelf.m']);
        if ~isempty(strfind(h, char([13, 10])))
            h(h == char(13)) = [];
        end
        h = splittocellc(h, char(10));
        hhelp = find(~cellfun('isempty', regexp(h, '^\%')));
        hhelp((1+findfirst(diff(hhelp) > 1)):end) = [];
        h = gluetostringc(regexprep(h(hhelp), '^\%\s?', ''), char(10), true);
        if nargout == 0
            disp(h);
        else
            varargout{1} = h;
        end
        return;
    end
    
    % invalid call to help
    if ~iscell(S(2).subs) || numel(S(2).subs) ~= 1 || ...
       ~ischar(S(2).subs{1}) || isempty(S(2).subs{1}) || ...
       ~isfield(neos.meth, S(2).subs{1})
        warning('neuroelf:general:badArgument', ...
            'Invalid help topic requested or topic not found.');
    end
    h = neos.meth.(S(2).subs{1}){end};
    if nargout == 0
        disp(h);
    else
        varargout{1} = h;
    end

    % return
    return;
end

% invalid method?
if ~isfield(neos.meth, S(1).subs)
    error('neuroelf:internal:badMethod', 'Invalid method or method not found.');
end

% check inputs/outputs (number)
meth = neos.meth.(S(1).subs);
if nargout > meth{2}
    warning('neuroelf:general:tooManyOutputs', 'Too many outputs requested.');
end
if numel(S) > 1
    if ~iscell(S(2).subs)
        error('neuroelf:general:badArguments', 'Invalid arguments in call.');
    end
    vargin = S(2).subs(:)';
    nnargin = numel(vargin);
else
    vargin = {};
    nnargin = 0;
end
if nnargin > meth{4}
    error('neuroelf:general:tooManyInputs', 'Too many inputs provided.');
end

% special case, one output but no inputs
if numel(S) == 1 && nargout == 1

    % return function handle instead
    varargout{1} = meth{1};
    return;
end

% try/catch
try

    % depending on call
    if meth{2} == 0 || nargout == 0
        feval(meth{1}, vargin{:});
    elseif nargout <= 1
        varargout{1} = feval(meth{1}, vargin{:});
    else
        [varargout{1:min(meth{2}, nargout)}] = feval(meth{1}, vargin{:});
    end

% error occurred
catch ne_eo;

    % rethrow
    rethrow(ne_eo);
end
