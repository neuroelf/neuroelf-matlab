function dispstruct(s)
% dispstruct  - display struct contents
%
% FORMAT:       dispstruct(s)
%
% Input fields:
%
%       s           struct variable
%
% No output fields.

% Version:  v0.9c
% Build:    11060413
% Date:     May-24 2011, 12:31 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Jochen Weber
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
   ~isstruct(s)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% get fieldnames and size
f = fieldnames(s);
z = size(s);

% no fields
if numel(f) == 0

    % simply print size
    zs = sprintf('%dx', z);
    disp([zs(1:end-1) ' struct array with no fields.']);

% fields
else

    % size not 1x1
    if ~isequal(z, [1, 1])

        % print size and fieldnames
        zs = sprintf('%dx', z);
        disp([zs(1:end-1) 'struct array with fields:']);
        for c = 1:numel(f)
            disp(['    ' f{c}]);
        end

    % size of struct is 1x1
    else

        % get maximum field length and add 4
        m = sprintf('%%%ds: %%s', size(char(f), 2) + 4);

        % iterate over fields
        for c = 1:numel(f)

            % get field content and size
            t = s.(f{c});
            l = class(t);
            z = size(t);
            zs = sprintf('%dx', z);
            zs = sprintf('[%s %s]', zs(1:end-1), l);

            % depending on type
            switch (lower(l))

                % char array
                case {'char'}

                    % not a (short) line
                    if ndims(t) ~= 2 || ...
                        size(t, 1) ~= 1 || ...
                        numel(t) > 48 || ...
                       ~all(t >= 32 & t <= 126)

                        % special case: empty string
                        if all(size(t) == 0)
                            t = '''''';

                        % otherwise
                        else

                            % replace with class
                            t = zs;
                        end

                    % otherwise, pack into quotes
                    else
                        t = sprintf('''%s''', t);
                    end

                % numeric array
                case {'double', 'single'}

                    % real values
                    if isreal(t)

                        % single number
                        if numel(t) == 1
                            t = sprintf('%g', t(1));

                        % short row of numbers
                        elseif ndims(t) == 2 && ...
                            size(t, 1) == 1 && ...
                            numel(t) <= 8 && ...
                           ~isempty(t)
                            t = sprintf('%g ', t);

                        % empty array
                        elseif all(z == 0)
                            t = '[]';

                        % otherwise
                        else
                            t = zs;
                        end

                    % complex content
                    else

                        % single number
                        if numel(t) == 1
                            t = sprintf('%g + %gi', real(t), complex(t));

                        % otherwise
                        else
                            zs = sprintf('%dx', z);
                            t = sprintf('%s %s (complex)', zs(1:end-1), l);
                        end
                    end

                % boolean values
                case {'logical'}

                    % one value
                    if numel(t) == 1

                        % value
                        if t
                            t = '[true]';
                        else
                            t = '[false]';
                        end

                    % multiple values
                    else
                        t = zs;
                    end

                % numeric value
                case {'int64', 'int32', 'int16', 'int8', ...
                      'uint64', 'uint32', 'uint16', 'uint8'}

                    % one value
                    if numel(t) == 1

                        % print value
                        t = sprintf('%s(%d)', l, t(1));

                    % otherwise
                    else
                        t = zs;
                    end

                % cell array
                case {'cell'}
                    t = zs;

                % unknown type
                otherwise

                    % print class string
                    t = zs;
            end

            % output
            disp(sprintf(m, f{c}, t));
        end
    end
end
