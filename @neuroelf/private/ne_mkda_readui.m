% PUBLIC FUNCTION ne_mkda_readui: test UI and read to analysis
function varargout = ne_mkda_readui(varargin)

% Version:  v0.9c
% Build:    12110111
% Date:     Nov-08 2011, 1:08 PM EST
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

% global variable
global ne_gcfg;
hFig = ne_gcfg.h.MKDA.MKDAFig;
ch = ne_gcfg.h.MKDA.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% action string
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3}) || ...
   ~any(strcmp(varargin{3}(:)', {'mask', 'nulspat', 'nulunit', 'scgauss', 'scindic'}))
    return;
end

% depending on action
switch (varargin{3}(:)')

    % mask selection
    case {'mask'}

        % only update on second index
        if ch.Mask.Value ~= 2
            return;
        end

        % set value back
        ch.Mask.Value = 1;

        % get current selection
        mskstring = ch.Mask.String;
        if ~iscell(mskstring)
            mskstring = cellstr(mskstring);
        end
        mskstring = mskstring{1};

        % get path
        colinpath = neuroelf_path('colin');
        if any(mskstring) == '/'
            mskpath = fileparts(mskstring);
        else
            mskpath = colinpath;
        end
        opwd = pwd;
        if exist(mskpath, 'dir') == 7
            cd(mskpath);
        end

        % select new mask
        [nmskfile, nmskpath] = uigetfile( ...
            {'*.vmr;*.msk;*.hdr;*.nii', ...
             'MKDA-compatible masking files (*.vmr, *.msk, *.hdr, *.nii)'}, ...
            'Please select a different mask file...');

        % switch back
        cd(opwd);

        % output
        if isequal(nmskfile, 0) || ...
            isequal(nmskpath, 0) || ...
            isempty(nmskfile) || ...
            isempty(regexpi(nmskfile, '^.*\.(hdr|msk|nii|vmr)$'))
            nmskfile = '<no mask selected>';
        else

            % requires path
            nmskpath = strrep(nmskpath, '\', '/');
            if ~isempty(nmskpath) && ...
                nmskpath(end) == '/'
                nmskpath(end) = '';
            end
            if ~strcmpi(nmskpath, colinpath)
                nmskfile = [nmskpath '/' nmskfile];
            end
        end

        % update mask file
        mskstring = {nmskfile; '<select HDR/MSK/VMR file>'};
        ch.Mask.String = mskstring;

        % list points
        ne_mkda_listpoints;

    % spatial null distribution
    case {'nulspat'}

        % update radio button group
        hFig.RadioGroupSetOne('SUNull', 1);

        % test for differential contrast
        if any(ch.Contrast.String == '>')

            % warning
            uiwait(warndlg( ...
                'Differential contrasts might not be valid with a spatial null!', ...
                'NeuroElf - warning', 'modal'));
        end

    % unit-based null
    case {'nulunit'}

        % update radio button group
        hFig.RadioGroupSetOne('SUNull', 2);

        % test for differential contrast
        if ~any(ch.Contrast.String == '>')

            % warning
            uiwait(warndlg( ...
                ['Unit label null hypotheses are most appropriate for differential ' ...
                 'contrasts (testing function not spatial specificity)!'], ...
                'NeuroElf - warning', 'modal'));
        end

    % sphere type
    case {'scgauss'}
        hFig.RadioGroupSetOne('Scale', 2);
        hFig.SetGroupEnabled('Gauss', 'on');
    case {'scindic'}
        hFig.RadioGroupSetOne('Scale', 1);
        hFig.SetGroupEnabled('Gauss', 'off');
end

% update analysis
ne_mkda_updana;
