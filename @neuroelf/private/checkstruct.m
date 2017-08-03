function cstr = checkstruct(cstr, varargin)
% checkstruct  - check an 1x1 struct for validity
%
% FORMAT:       checked = checkstruct(tocheck, options [, remunknown])
%
% Input fields:
%
%       tocheck     1x1 struct to check
%       options     Nx4 cell array with options, whereas for each row
%                   has the form or {field, classname, condition, default}, and
%                   the first cell contains the field name
%                   the second cell is one out of
%                   {'char', 'numeric', 'struct', 'cell'}
%                   the third cell is one out of
%                   {'nonempty', 'noinfnan', 'label', 'expression', 'deblank'}
%                   and the last cell is the default for when either the
%                   type or condition doesn't match
%       remunknown  if given removes all unnamed fields from tocheck
%
% Output fields:
%
%        checked    struct after type checking

% Version:  v1.1
% Build:    16042016
% Date:     Apr-20 2016, 4:25 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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

% basic argument check
if nargin < 2 || ...
  ~iscell(varargin{1}) || ...
   isempty(varargin{1}) || ...
   size(varargin{1}, 2) ~= 4
    error( ...
        'neuroelf:BadArgument', ...
        'Checkstruct needs a Mx4 options argument.' ...
    );
end
opt = varargin{1};

% build struct from non structs
if ~isstruct(cstr) || ...
    length(cstr) ~= 1
    cstr = struct;
end

% remove unknown fields
if nargin > 2
    remunknown = true;
    flist = struct;
else
    remunknown = false;
end

% grand try block for bad sub-types
try

    % iterate over fields
    for fc = 1:size(opt, 1)
        tname = opt{fc, 1}(:)';

        % keep track of fieldnames
        if remunknown
            flist.(tname) = true;
        end

        % check field existance
        if ~isfield(cstr, tname)
            cstr.(tname) = opt{fc, 4};
            continue;
        end

        % get current content
        tcont = cstr.(tname);

        % check field class
        tclass = lower(opt{fc, 2}(:)');
        if isempty(tclass)
            continue;
        end
        tclcmp = lower(class(tcont));
        if isempty(regexpi(tclcmp, tclass))
            cstr.(tname) = opt{fc, 4};
            continue;
        end

        % options
        topt = opt{fc, 3}(:)';

        % character option
        if ischar(topt)

            % what option
            switch (lower(opt{fc, 3}(:)'))

                % deblanking
                case {'deblank'}
                    if ischar(tcont)
                        cstr.(tname) = deblank(tcont);
                    end

                % expression
                case {'expression'}
                    if ischar(tcont) && ~isempty(checksyntax(tcont(:)'))
                        cstr.(tname) = opt{fc, 4};
                    elseif ischar(tcont)
                        cstr.(tname) = deblank(tcont(:)');
                    end

                % label
                case {'label'}
                    if ~isrealvarname(deblank(tcont(:)'))
                        cstr.(tname) = opt{fc, 4};
                    else
                        cstr.(tname) = deblank(tcont(:)');
                    end

                % nonempty
                case {'nonempty'}
                    if isempty(cstr.(tname))
                        cstr.(tname) = opt{fc, 4};
                    end

                % noinfnan
                case {'noinfnan'}
                    if isnumeric(tcont) && ...
                       (any(isinf(tcont)) || any(isnan(tcont)))
                        cstr.(tname) = opt{fc, 4};
                    end

                % singular content
                case {'singular'}
                    if numel(cstr.(tname)) ~= 1
                        cstr.(tname) = opt{fc, 4};
                    end
            end

        % cell option
        elseif iscell(topt)

            % test content
            cmatch = false;
            tsize  = numel(tcont);
            for cc = 1:numel(topt)
                if numel(topt{cc}) == tsize && ...
                    all(topt{cc} == tcont)
                    cmatch = true;
                    break;
                end
            end

            % match not found
            if ~cmatch
                cstr.(tname) = opt{fc, 4};
            end
        end
    end

catch ne_eo;
    error( ...
        'neuroelf:BadArgument', ...
        'A suboption of options seems invalid (%s).', ...
        ne_eo.message ...
    );
end
