function [varargout] = ne_vtc_filter(varargin)
% ne_vtc_filter  - apply filtering to VTC content
%
% FORMAT:       [vtc = ] ne_vtc_filter(SRC, EVT, type [, opts])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       type        type of filtering, 'temp' or 'spat'
%       opts        optional settings
%        .tcutoff   temporal cut-off, 1/frequency (sec, e.g. 128)
%        .ttype     temporal filter type, see neuroelf::tempfilter
%
% Output fields:
%
%       vtc         filtered VTC (in place!)
%
% Notes: the VTC content will be filtered in place (transio loaded);
%        requires the current SliceVar object to be of type VTC!

% Version:  v1.0
% Build:    15031017
% Date:     Mar-10 2015, 5:27 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2015, Jochen Weber
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

% preset output
varargout = cell(1, nargout);

% global variable
global ne_gcfg;
vtc = ne_gcfg.fcfg.SliceVar;
if numel(vtc) ~= 1 || ...
   ~isxff(vtc, 'vtc')
    return;
end

% create VMP
if nargout > 0
    varargout{1} = vtc;
end

% nothing to do
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3}) || ...
   ~any(strcmpi(varargin{3}(:)', {'spat', 'temp'}))
    return;
end

% options
if nargin < 4 || ...
   ~isstruct(varargin{4}) || ...
    numel(varargin{4}) ~= 1
    opts = struct;
else
    opts = varargin{4};
end

% spatial
if lower(varargin{3}(1)) == 's'
    
% temporal
else
    
    % options given
    if nargin > 3 && ...
        ischar(varargin{4}) && ...
       ~isempty(varargin{4}) && ...
        any(strcmpi(varargin{4}(:)', {'dct', 'fourier', 'poly'}))
        opts.ttype = lower(varargin{4}(:)');
        if nargin > 4 && ...
            ischar(varargin{5}) && ...
           ~isempty(varargin{5}) && ...
           ~isempty(regexpi(ddeblank(varargin{5}(:)'), '^\d+(\.\d+)?$'))
            opts.tcutoff = str2double(ddeblank(varargin{5}(:)'));
        end
    end

    % not all options given/valid
    if ~isfield(opts, 'ttype') || ...
       ~ischar(opts.ttype) || ...
        isempty(opts.ttype) || ...
       ~any(strcmpi(opts.ttype, {'dct', 'fourier', 'poly'})) || ...
       ~isfield(opts, 'tcutoff') || ...
       ~isa(opts.tcutoff, 'double') || ...
        numel(opts.tcutoff) ~= 1 || ...
        isinf(opts.tcutoff) || ...
        isnan(opts.tcutoff) || ...
        opts.tcutoff <= (0.004 * vtc.TR)
        
        % get options
        topts = {'  dct', '  128'};
        if isfield(opts, 'ttype') && ...
            ischar(opts.ttype) && ...
            any(strcmpi(opts.ttype, {'fourier', 'poly'}))
            topts{1} = ['  ' lower(opts.ttype)];
        end
        if isfield(opts, 'tcutoff') && ...
            isa(opts.tcutoff, 'double') && ...
            numel(opts.tcutoff) == 1 && ...
           ~isinf(opts.tcutoff) && ...
           ~isnan(opts.tcutoff) && ...
            opts.tcutoff >= (0.004 * vtc.TR)
            topts{2} = sprintf('  %d', 0.1 * round(0.04 * vtc.TR));
        end
        topts = inputdlg({'Filtering basis set (dct, fourier, poly):', ...
            'Filtering cutoff (1/Hz, seconds):'}, 'NeuroElf - user input', ...
            1, topts);
        if numel(topts) ~= 2 || ...
            isempty(topts{1}) || ...
           ~any(strcmpi(ddeblank(topts{1}), {'dct', 'fourier', 'poly'})) || ...
            isempty(topts{2}) || ...
            isempty(regexpi(ddeblank(topts{2}), '^\d+(\.\d+)?$'))
            return;
        end
        opts.ttype = ddeblank(topts{1});
        opts.tcutoff = str2double(ddeblank(topts{2}));
    end
    
    % nothing to be done
    if opts.tcutoff >= (0.001 * vtc.TR * vtc.NrOfVolumes)
        return;
    end
    
    % implement filter
    opts.temp = true;
    if opts.ttype(1) == 'd'
        opts.tempdct = opts.tcutoff;
    elseif opts.ttype(1) == 'f'
        opts.tempsc = floor(sqrt(eps) + (0.001 * vtc.TR * vtc.NrOfVolumes) / opts.tcutoff);
    else
        opts.temppoly = max(1, floor(sqrt(eps) + (0.002 * vtc.TR * vtc.NrOfVolumes) / opts.tcutoff));
    end
    
    % UI and filter
    mfp = ne_gcfg.h.MainFig.Pointer;
    ne_gcfg.h.MainFig.Pointer = 'watch';
    drawnow;
    vtc.Filter(opts);
    ne_gcfg.h.MainFig.Pointer = mfp;
    ne_setslicepos;
end
