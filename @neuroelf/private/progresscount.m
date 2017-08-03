function [varargout] = progresscount(varargin)
% progresscount  - display progress on a step based counter
%
% FORMAT:       progresscount(totalsteps, dotsize, percentsize, character)
%
% Input fields:
%
%       totalsteps  number of items/steps being processed
%       dotsize     percents needed for a "dot" to print out
%       percentsize percents needed for printing a percentage
%       character   character to use as a "dot" character
%
% Example:
%
%       progresscount(40000, 2.5, 10, '.');
%
% initializes the routine to print out a dot every 2.5 percent,
% and a percentage every 10 percent. progress character is '.'
%
%
% FORMAT:       progresscount(currentstep)
%
% Input fields:
%
%       currentstep number of current step being processed
%
% Example:
%
%       progresscount(14200);
%
% will print all necessary characters, percentages for this step
%
%
% FORMAT:       whereami = progresscount;
%
% returns the current progress state as a real number [0.0;1.0]

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

% persistent progress storage
persistent iv_progress;
if isempty(iv_progress)
    iv_progress = struct('init', 0);
end

% what to do
switch nargin

    % no input
    case {0}

        % output
        if nargout > 0

            % whereami only if init
            if iv_progress.init == 1
                varargout{1} = iv_progress.count / iv_progress.total;

            % otherwise 0
            else
                varargout{1} = 0;
            end
        end

    % one input
    case 1

        % only valid if init
        if iv_progress.init ~= 1
            if nargout > 0
                varargout{1} = 0;
            end
            return;
        end

        % only for numeric values as well
        if ~isnumeric(varargin{1}) || ...
            isempty(varargin{1}) || ...
            isinf(varargin{1}(1)) || ...
            isnan(varargin{1}(1)) || ...
            varargin{1}(1) < 0
            return;
        end

        % get counter
        count = varargin{1}(1);

        % set internal value
        iv_progress.count = min(count, iv_progress.total);

        % while below next dot
        while (100 * (count / iv_progress.total) > iv_progress.dsnext)

            % set next
            iv_progress.dsnext = iv_progress.dsnext + iv_progress.step;

            % give percentage
            if iv_progress.dsnext >= iv_progress.psnext

                % show percentage
                fprintf(1, ' %3.0f%% ', iv_progress.psnext);

                % update percentage limit
                iv_progress.psnext = iv_progress.psnext + iv_progress.pstep;
                continue;
            end

            % otherwise print char
            fprintf(1, iv_progress.char);
        end

        % give percentage?
        if nargout == 1
            varargout{1} = count / iv_progress.total;
        end

    % init syntax
    otherwise

        % only for 4 arguments
        if nargin < 4
            return;
        end

        % set arguments
        iv_progress.total = max(varargin{1}, 1);
        iv_progress.step  = varargin{2};
        iv_progress.pstep = varargin{3};
        iv_progress.char  = varargin{4};
        iv_progress.dsnext= 0;
        iv_progress.psnext= 0;
        iv_progress.count = 0;
        iv_progress.init = 1;
end
