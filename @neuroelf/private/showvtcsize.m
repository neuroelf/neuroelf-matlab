function [varargout] = showvtcsize(varargin)
% showvtcsize  - displays the sizes of VTCs found on the given path
%
% FORMAT:       [s] = showvtcsize(folder, [pattern, ...])
%
% Input fields:
%
%       folder      folder where to look for VTCs
%       pattern     pattern
%       ...         options passed to findfiles
%
% Output fields:
%
%       s           1xN struct with fields
%        .filename  VTC filename
%        .protocol  linked protocol
%        .res       VTC resolution
%        .size      VTCData size
%        .xyzstart  start coordinate

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-08 2011, 9:16 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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

% initialize output first
s = emptystruct({'filename', 'protocol', 'res', 'size', 'xyzstart'}, [1, 0]);

% check arguments
if nargin < 1 || ...
   ~ischar(varargin{1}) || ...
    isempty(varargin{1})
    fld = pwd;
    v = findfiles(fld, '*.vtc', 'relative=');
elseif exist(varargin{1}(:)', 'dir') == 7
    fld = varargin{1}(:)';
    if nargin == 1 || ...
       ~ischar(varargin{2}) || ...
        isempty(varargin{2})
        try
            v = findfiles(fld, '*.vtc', 'relative=', varargin{2:end});
        catch ne_eo;
            rethrow(ne_eo);
        end
    else
        if any(varargin{2}(:) == '*')
            try
                v = findfiles(fld, varargin{2}, 'relative=', varargin{3:nargin});
            catch ne_eo;
                rethrow(ne_eo);
            end
        else
            try
                v = findfiles(fld, 'relative=', varargin{2:nargin});
            catch ne_eo;
                rethrow(ne_eo);
            end
        end
    end
elseif any(varargin{1}(:) == '*')
    fld = pwd;
    try
        v = findfiles(fld, 'relative=', varargin{:});
    catch ne_eo;
        rethrow(ne_eo);
    end
else
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid pattern or directory name given.' ...
    );
end

% check output
for vc = numel(v):-1:1
    if isempty(regexpi(v{vc}, '\.vtc$'))
        v(vc) = [];
    end
end
if isempty(v)
    return;
end

% iterate over files
s(numel(v)).filename = v{end};
for vc = 1:numel(s)

    % set filename
    s(vc).filename = v{vc};

    % load VTC
    try
        vtc = xff([fld '/' v{vc}], 'h');
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        s(vc).filename = '';
        continue;
    end

    % get settings
    s(vc).protocol = vtc.NameOfLinkedPRT;
    s(vc).res = vtc.Resolution;
    s(vc).size = size(vtc.VTCData);
    s(vc).xyzstart = [vtc.XStart, vtc.YStart, vtc.ZStart];
end

% clean output list and get max name sizes
mxf = 8;
mxp = 8;
for vc = numel(s):-1:1
    if isempty(s(vc).filename)
        s(vc) = [];
        continue;
    end
    mxf = max(mxf, numel(s(vc).filename));
    mxp = max(mxp, numel(s(vc).protocol));
end

% show sizes or return with struct
if nargout < 1
    mxf = min(mxf, 56);
    mxp = min(mxp, 56);
    mxf = sprintf(' %%-%ds | %%-%ds | ', mxf, mxp);
    sout = sprintf( ...
        [mxf 'Volumes | Dimensions  | Res | Start'], 'Filename', 'Protocol');
    disp(sout);
    sout(:) = '-';
    disp(sout);
    mxf = [mxf '%4d    | %3d,%3d,%3d | %2d  | %3d,%3d,%3d'];
    for vc = 1:numel(s)
        fn = s(vc).filename;
        if numel(fn) > 56
            fn = [fn(1:24) '...' fn(end-24:end)];
        end
        pt = s(vc).protocol;
        if numel(pt) > 56
            pt = [pt(1:24) '...' pt(end-24:end)];
        end
        disp(sprintf(mxf, fn, pt, s(vc).size(1), ...
            s(vc).size(2), s(vc).size(3), s(vc).size(4), s(vc).res, ...
            s(vc).xyzstart(1), s(vc).xyzstart(2), s(vc).xyzstart(3)));
    end
else
    varargout = cell(1, nargout);
    varargout{1} = s;
end
