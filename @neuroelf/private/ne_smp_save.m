% FUNCTION ne_smp_save: save currently selected SurfStatsVar
function varargout = ne_smp_save(varargin)

% Version:  v0.9d
% Build:    14071116
% Date:     Jul-11 2014, 4:39 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, Jochen Weber
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
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% get/check current SMP
smp = cc.SurfStatsVar;
if numel(smp) ~= 1 || ...
   ~isxff(smp, 'smp')
    return;
end

% no further arguments
e = [];
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    numel(varargin{3}) < 4

    % echo
    if ne_gcfg.c.echo
        ne_echo('smp', 'Save');
    end

    % try saving to disk
    try
        if ~isempty(smp.FilenameOnDisk)
            smp.Save;
        else
            smp.SaveAs;
        end
    catch ne_eo;
        e = ne_eo;
    end

% filename given
elseif strcmpi(varargin{3}(end-3:end), '.smp')

    % echo
    if ne_gcfg.c.echo
        ne_echo('smp', sprintf('SaveAs(%s)', varargin{3}(:)'))
    end

    % save
    try
        smp.SaveAs(varargin{3}(:)');
    catch ne_eo;
        e = ne_eo;
    end

% general "saveas" call
elseif strcmpi(varargin{3}(:)', 'saveas')

    % echo
    if ne_gcfg.c.echo
        ne_echo('smp', 'SaveAs');
    end

    % try saving as
    try
        smp.SaveAs;
    catch ne_eo;
        e = ne_eo;
    end
end

% error handing
if ~isempty(e)
    uiwait(warndlg(e.message, 'NeuroElf - warning', 'modal'));

% otherwise
else

    % update filename
    sidx = ch.SurfStatsVar.Value;
    svars = ch.SurfStatsVar.UserData;
    [fp, fn, fe] = fileparts(smp.FilenameOnDisk);
    fn = [fn, fe];
    if isempty(fe)
        return;
    end
    if sidx <= size(svars, 1) && ...
        isxff(svars{sidx, 4}, 'smp') && ...
        smp == svars{sidx, 4}
        svarnames = ch.SurfStatsVar.String;
        if ischar(svarnames)
            svarnames = cellstr(svarnames);
        end
        svarnames{sidx} = fn;
        ch.SurfStatsVar.String = svarnames;
        ch.SurfStatsVar.UserData{sidx, 1} = fn;
    end
end
