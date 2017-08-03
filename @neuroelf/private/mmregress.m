function [ixx, varargout] = mmregress(x, varargin)
% mmregress  - regress multiple models / multiple data
%
% FORMAT:       [ixx, b, ...] = mmregress(x, y [, ...])
%
% Input fields:
%
%       x           NxR or NxRxM model(s)
%       y, ...      Nx1 or NxM (or Nx1xM) data
%
% Output fields:
%
%       ixx         inverse of the x' * x matrix/ces
%       b, ...      regression beta estimates
%
% Note: please use this function only for a limited set of regressors,
%       as the precision will severely suffer, also if regressors are
%       highly correlated, given the use of the pseudo-inverse function

% Version:  v1.0
% Build:    14091915
% Date:     Sep-19 2014, 3:48 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2013, 2014, Jochen Weber
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

% argument check (do only the most basic things)
if nargin < 2 || ...
   ~isa(x, 'double') || ...
    ndims(x) > 3 || ...
   (~isa(varargin{1}, 'double') && ...
    ~isa(varargin{1}, 'single')) || ...
    size(x, 1) ~= size(varargin{1}, 1)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% set output size
nvi = nargin - 1;
varargout = cell(1, nvi);

% with error handling...
try

    % get the model size
    sx = size(x);

    % if the model is a Xx1 vector
    if sx(2) == 1

        % the inverse of the model is simply 1 / sum(x .* x)
        ixx = 1 ./ sum(x .* x);

        % adapt the model for multiplication with the data
        if numel(ixx) == 1
            x = reshape(ixx .* x, 1, sx(1));
        else
            x = reshape(ixx(ones(sx(1), 1), :, :) .* x, sx([2, 1, 3]));
        end

        % iterate over additional inputs
        for ac = 1:nvi

            % get data and size
            d = varargin{ac};
            sd = size(d);

            % model is the same for all data
            if numel(ixx) == 1

                % regress the single model directly using matlab's mtimes
                p = reshape(x * reshape(d, sd(1), prod(sd(2:end))), [sx(2), sd(2:end)]);

            % model is multi-dim, and data also
            elseif numel(d) > sx(1)

                % attempt using transmul (for models and then for data)
                p = reshape(transmul(x, ...
                    reshape(d, [sd(1), 1, prod(sd(2:end))])), [sx(2), sd(2:end)]);

            % model is multi-dim, but only one set of data
            else

                % replicate data when multiplying with models
                p = reshape(transmul(x, ...
                    repmat(d, [1, 1, sx(3)])), [sx(2), 1, sx(3:end)]);
            end

            % add to output
            varargout{ac} = p;
        end
        return;
    end

    % inverse of model(s)
    ixx = invnd(transmul(x));

    % repeat for each additional input
    for ac = 1:nvi

        % get data
        d = varargin{ac};

        % size of data
        sd = size(d);

        % model is the same for all data
        if numel(sx) == 2

            % regress the single model directly using matlab's mtimes
            p = reshape((ixx * x') * ...
                reshape(d, sd(1), prod(sd(2:end))), [sx(2), sd(2:end)]);

        % model is multi-dim, and data also
        elseif numel(d) > sx(1)

            % attempt using transmul (for models and then for data)
            p = reshape(transmul(transmul(ixx, x, 2), ...
                reshape(d, [sd(1), 1, prod(sd(2:end))])), [sx(2), sd(2:end)]);

        % model is multi-dim, but only one set of data
        else

            % replicate data when multiplying with models
            p = reshape(transmul(transmul(ixx, x, 2), ...
                repmat(d, [1, 1, sx(3)])), [sx(2), 1, sx(3:end)]);
        end

        % add to output
        varargout{ac} = p;
    end

% pass on errors
catch ne_eo;
    rethrow(ne_eo);
end
