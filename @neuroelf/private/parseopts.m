function opts = parseopts(opts, varargin)
% parseopts  - parse optional input arguments with defaults
%
% FORMAT:       opts = parseopts(opts, ...)
%
% Input fields:
%
%       opts        1x1 struct with defaults for options
%       ...         string/value pairs, 1x1 structs or opt=value strings
%
% Output fields:
%
%       opts        1x1 struct with filled options
%
% Note: input arguments are processed in given order;
%       argument names in ... are case *insensitive*

% Version:  v0.9a
% Build:    10051716
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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

% argument check
if nargin < 1 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1 || ...
    length(fieldnames(opts)) < 1
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing opts struct.' ...
    );
end

% get list of valid options and args
hasopts = fieldnames(opts);
valopts = lower(hasopts);
valvars = varargin;

% parse all other args
ac = 1;
while ac <= length(valvars)

    % get next argument
    nextarg = valvars{ac};
    ac = ac + 1;

    % reject empty arguments
    if isempty(nextarg)
        continue;
    end

    % type of argument
    switch (lower(class(nextarg)))

        % char arguments
        case {'char'}

            % linearize
            nextarg = nextarg(:)';

            % is a name=val option
            nameval = find(nextarg == '=');
            if ~isempty(nameval)
                valvars = [valvars(1:ac-1), ...
                    {nextarg(1:nameval(1)-1), nextarg(nameval(1)+1:end)}, ...
                    valvars(ac:end)];

            % option, value pair
            else

                % ac invalid already
                if ac > length(valvars)
                    error( ...
                        'neuroelf:TooFewArguments', ...
                        'Missing value for option ''%s''.', ...
                        nextarg ...
                    );
                end

                % get option value and increase ac
                nextval = valvars{ac};
                ac = ac + 1;

                % find correct option
                parnum = find(strcmpi(nextarg, valopts));
                if isempty(parnum)
                    warning( ...
                        'neuroelf:UnknownOption', ...
                        'Unknown option name: ''%s''.', ...
                        nextarg ...
                    );
                    continue;
                end

                % set option
                opts.(hasopts{parnum(1)}) = nextval;

            end

        % struct arguments
        case {'struct'}

            % only valid for 1x1 structs
            if numel(nextarg) ~= 1
                continue;
            end

            % put contents into valvars
            argf = fieldnames(nextarg);
            argc = {};
            for fc = 1:length(argf)
                argc = [argc, {argf{fc}, nextarg.(argf{fc})}];
            end
            valvars = [valvars(1:ac-1), argc(:)', valvars(ac:end)];

        % unsupported type
        otherwise
            warning( ...
                'neuroelf:BadArgument', ...
                'Unsupported input argument type.' ...
            );
            continue;
    end
end
