% FUNCTION ne_updatename: update name
function varargout = ne_updatename(varargin)

% Version:  v1.0
% Build:    16010921
% Date:     Jan-09 2016, 9:24 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2015, 2016, Jochen Weber
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

% global access to configuration/handles
global ne_gcfg;
cc = ne_gcfg.c;
ch = ne_gcfg.h;

% prepare output
varargout = cell(1, nargout);

% add progress
if cc.progress{1}
    if ~isempty(cc.progress{3})
        task = cc.progress{3};
        if numel(task) > 24
            spacepos = findfirst(task == ' ', 24 - numel(task));
            if isempty(spacepos) || ...
                spacepos == 1
                task = [task(1:21) '...'];
            else
                task = task(1:(spacepos-1));
            end
        end
    else
        task = 'task';
    end
    progress = ddeblank(sprintf('%5.1f%%', 100 * cc.progress{2}));
    progress = sprintf('[%s %s] - ', progress, task);
else
    progress = '';
end

% page determines title row
if ne_gcfg.fcfg.page == 3
    titlerow = 2;
else
    titlerow = 1;
end
titles = cc.title(titlerow, :);

% compose from title particles
if ~isempty(titles{2})
    if ~isempty(titles{3})
        title = sprintf('NeuroElf GUI - %s%s - %s :: %s', progress, titles{1:3});
    else
        title = sprintf('NeuroElf GUI - %s%s - %s', progress, titles{1:2});
    end
else
    title = sprintf('NeuroElf GUI - %s%s', progress, titles{1});
end

% set name and update
ch.MainFig.Name = title;
if nargin < 3 || ...
   ~islogical(varargin{3}) || ...
    numel(varargin{3}) ~= 1 || ...
    varargin{3}
    drawnow;
end
