function xo = ntt_CleanupStimChannel(xo, chnr, opts)
% NTT::CleanupStimChannel  - perform stim channel cleanup
%
% FORMAT:       [ntt = ] ntt.CleanupStimChannel(chnr, opts);
%
% Input fields:
%
%       chnr        number of the channel containing the stim info
%       opts        settings
%        .cmd       mandatory 1xC cell with steps to perform, selection of
%                   'renum' - renumber trains of values in channel
%        .renum     Rx2 cell array with renumbering trains (mandatory)
%        .rminint   mininum change interval (default: 2 samples) for renum
%
% Output fields:
%
%       ntt         altered object
%
% Using: findfirst.

% Version:  v1.1
% Build:    16021215
% Date:     Feb-12 2016, 3:34 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2014, 2016, Jochen Weber
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

% neuroelf library
global ne_methods;
findfirst = ne_methods.findfirst;

% argument check
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, {'acq', 'ntt'}) || ...
    numel(chnr) ~= 1 || ~isa(chnr, 'double') || isinf(chnr) || isnan(chnr) || ...
    chnr < 1 || chnr ~= fix(chnr) || numel(opts) ~= 1 || ...
   ~isstruct(opts) || ~isfield(opts, 'cmd') || isempty(opts.cmd)
    error('neuroelf:xff:badArgument', 'Bad or missing argument for call.');
end
if ischar(opts.cmd)
    opts.cmd = {opts.cmd(:)'};
elseif ~iscell(opts.cmd)
    error('neuroelf:xff:badArgument', 'Bad of missing .cmd option.');
end

% get content and filetype
bc = xo.C;
ftype = lower(xo.S.Extensions{1});

% get the requested data
switch (ftype)

    % from ACQ objects
    case 'acq'

        % later

    % from NTT objects
    case 'ntt'

        % in the .Data field
        data = bc.Data;

        % check channel number for bounds
        if chnr > size(data, 2)
            error('neuroelf:xff:badArgument', 'Channel number out of bounds.');
        end

        % and get only this data
        data = data(:, chnr);
end

% iterate over steps
for cmc = 1:numel(opts.cmd)

    % not a valid command
    if ~ischar(opts.cmd{cmc}) || ~any(strcmpi(opts.cmd{cmc}(:)', {'renum'}))
        warning('neuroelf:xff:unsupported', ...
            'Unsupported command at position %d in list.', cmc);
        continue;
    end

    % what command
    switch (lower(opts.cmd{cmc}(:)'))

        % renumbering
        case 'renum'

            % parse partial options
            if ~isfield(opts, 'renum') || ~iscell(opts.renum) || ...
                ndims(opts.renum) > 2 || isempty(opts.renum) || size(opts.renum, 2) ~= 2
                error('neuroelf:xff:badArgument', ...
                    'Bad or missing .renum option for command ''renum''.');
            end
            renum = opts.renum;
            relen = zeros(size(renum, 1), 1);
            for rc = 1:size(renum, 1)
                if ~isa(renum{rc, 1}, 'double') || ~isa(renum{rc, 2}, 'double') || ...
                    numel(renum{rc, 1}) ~= numel(renum{rc, 1})
                    error('neuroelf:xff:badArgument', ...
                        'Stim trains must be double vectors of same length.');
                end
                renum{rc, 1} = renum{rc, 1}(:);
                renum{rc, 2} = renum{rc, 2}(:);
                relen(rc) = numel(renum{rc, 1});
            end
            [relen, rells] = sort(relen, 'descend');
            renum = renum(rells, :);
            if ~isfield(opts, 'rminint') || ~isa(opts.rminint, 'double') || numel(opts.rminint) ~= 1 || ...
                isinf(opts.rminint) || isnan(opts.rminint) || rminint < 1
                rminint = 2;
            else
                rminint = ceil(opts.rminint);
            end

            % get the positions of change in channel
            stchange = 1 + find(diff(data) ~= 0);
            stchange(end + 1) = numel(data) + 1;

            % remove those that are too short to be considered real
            stchange(find(diff(stchange) < rminint)) = [];

            % and get values
            stimdata = data(stchange(1:end-1));

            % now iterate over renum trains
            rmatch = zeros(numel(stimdata), size(renum, 1));
            for rc = 1:size(renum, 1)

                % locate places where this train matches
                tm = true(numel(stimdata) + 1 - numel(renum{rc, 1}), 1);
                for cc = 1:numel(renum{rc, 1})
                    tm = tm & (stimdata(cc:end+cc-numel(renum{rc, 1})) == renum{rc, 1}(cc));
                end
                if ~any(tm)
                    warning('neuroelf:xff:notFound', ...
                        'Train %d ( %s) not found.', rc, sprintf('%d ', renum{rc, 1}));
                    continue;
                end
                rmatch(find(tm), rc) = rc;
            end

            % iterate over stim data
            stc = 1;
            while stc <= numel(stimdata)

                % no train matched at this position
                tmatch = findfirst(rmatch(stc, :) > 0);
                if isempty(tmatch)

                    % next position
                    stc = stc + 1;
                    continue;
                end

                % assume the first train that matches is *the one*
                for rc = 1:numel(renum{tmatch, 1})

                    % replace values in data!
                    data(stchange(stc):stchange(stc+1)-1) = renum{tmatch, 2}(rc);
                    stc = stc + 1;
                    if stc > numel(stimdata)
                        break;
                    end
                end
            end
    end
end

% set data back
switch (lower(ftype))

    % for ACQ objects
    case 'acq'

        % later

    % for NTT objects
    case 'ntt'

        % place back into .Data field
        bc.Data(:, chnr) = data;
end

% and set back to storage
xo.C = bc;
