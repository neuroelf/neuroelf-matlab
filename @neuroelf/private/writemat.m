function writemat(WMI_NWV_matfile, varargin)
% writemat  - write MAT file with specific variables and content
%
% FORMAT:       writemat(matfile [, mfver], vname, vcont, ...)
%
% Input fields:
%
%       matfile     filename for target MAT file
%       mfver       valid version argument for save (e.g. '-v6')
%       vname/vcont pairs of names and content (as in struct(...))
%
% No output fields.
%
% Note: alternatively, vname and vcont can be cell arrays of the same
%       size, in which case mfver MUST be given however.

% Version:  v0.9d
% Build:    14061709
% Date:     Jun-17 2014, 9:45 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, Jochen Weber
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
if nargin < 3 || ...
   ~ischar(WMI_NWV_matfile) || ...
    isempty(WMI_NWV_matfile) || ...
   ~ischar(varargin{1}) || ...
    isempty(varargin{1})
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% internal variables begin with WMI_NWV_
if any(strcmpi(varargin{1}(:)', {'-v4', '-v6', '-v7', '-v7.3'}))
    WMI_NWV_mlversion = {lower(varargin{1}(:)')};
    WMI_NWV_firstarg = 2;
else
    WMI_NWV_mlversion = {};
    WMI_NWV_firstarg = 1;
end

% special case, input is given as 2 cell arrays
WMI_NWV_mmversion = mainver;
if nargin == 4 && ...
    iscell(varargin{2}) && ...
    iscell(varargin{3}) && ...
    isequal(size(varargin{2}), size(varargin{3})) && ...
    all(cellfun(@ischar, varargin{2}(:))) && ...
   ~any(cellfun('isempty', varargin{2}(:)))

    % assign in current workspace
    try
        for WMI_NWV_counter = 1:numel(varargin{2})
            eval(sprintf('%s=varargin{3}{WMI_NWV_counter};', ...
                varargin{2}{WMI_NWV_counter}(:)'));
        end
    catch ne_eo;
        rethrow(ne_eo);
    end

    % call save
    if WMI_NWV_mmversion > 6
        save(WMI_NWV_matfile(:)', varargin{2}{:}, WMI_NWV_mlversion{:}, '-mat');
    elseif WMI_NWV_mmversion > 4
        save(WMI_NWV_matfile(:)', varargin{2}{:}, '-mat');
    else
        if isempty(WMI_NWV_mlversion)
            WMI_NWV_mlversion{1} = '-v6';
        end
        save(WMI_NWV_mlversion{:}, WMI_NWV_matfile(:)', varargin{2}{:});
    end
else

    % assign variables to current workspace
    WMI_NWV_savevars = cell(1, ceil(0.5 * (nargin - WMI_NWV_firstarg)));
    try
        for WMI_NWV_counter = WMI_NWV_firstarg:2:(nargin-2)
            eval(sprintf('%s=varargin{WMI_NWV_counter+1};', varargin{WMI_NWV_counter}(:)'));
            WMI_NWV_savevars{ceil(0.5 * WMI_NWV_counter)} = varargin{WMI_NWV_counter}(:)';
        end
    catch ne_eo;
        rethrow(ne_eo);
    end

    % call save
    if WMI_NWV_mmversion > 6
        save(WMI_NWV_matfile(:)', WMI_NWV_savevars{:}, WMI_NWV_mlversion{:}, '-mat');
    elseif WMI_NWV_mmversion > 4
        save(WMI_NWV_matfile(:)', WMI_NWV_savevars{:}, '-mat');
    else
        if isempty(WMI_NWV_mlversion)
            WMI_NWV_mlversion{1} = '-v6';
        end
        save(WMI_NWV_mlversion{:}, WMI_NWV_matfile(:)', WMI_NWV_savevars{:});
    end
end
