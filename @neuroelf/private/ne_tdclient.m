function varargout = ne_tdclient(varargin)
% ne_tdclient  - TDclient interface for NeuroElf GUI
%
% FORMAT:       tdreply = ne_tdclient(SRC, EVT, [x, y, z], kind, cubesize)
% FORMAT:       ne_tdclient(SRC, EVT, 'visual')
%
% Input fields:
%
%       x, y, z     coordinates to be looked up
%       kind        string, {'TalLabel'}, 'Cube', 'NGM', 'RNGM'
%       cubesize    size of search cube or radius for NGM
%
% Output fields:
%
%       tdreply     reply of TD client, see help of tdclient function
%
% Example:
%
%   ne_tdclient(0, 0, 'visual');
%   reply = ne_tdclient(0, 0, [24, 0, -20], 'Cube')

% error handling
global ne_gcfg;

% simply pass on to tdclient
if nargin < 3
    tdclient('visual');
else
    try
        if nargout > 0
            [varargout{1:nargout}] = tdclient(varargin{3:end});
        else
            tdclient(varargin{3:end})
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
