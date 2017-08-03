function oppo = nepp_tosingle(ppo, opts)
% nepp_tosingle  - set imaging files to single precision
%
% FORMAT:       [oppo = ] nepp_tosingle(ppo [, opts])
%
% Input fields:
%
%       ppo         pre-processing object(s)
%       opts        optional settings
%        .prefix    if given and non-empty, create copy rather than in-place
%
% Output fields:
%
%       oppo        output pre-processing object(s) in same format
%
% Note: this function supports inputs of type
%       - char array or cell array of strings (filenames)
%       - DMR, FMR, HDR(NII), HEAD, and VTC objects

% Version:  v0.9d
% Build:    14081815
% Date:     Aug-18 2014, 3:10 PM EST
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

% argument check
if nargin < 1 || ...
   (~ischar(ppo) && ...
    ~iscell(ppo) && ...
    (numel(ppo) ~= 1 || ...
     ~isxff(ppo, {'dmr', 'fmr', 'hdr', 'head', 'vtc'})))
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid preprocessing object supplied.' ...
    );
end
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1 || ...
   ~isstruct(opts)
    opts = struct;
end
if ~isfield(opts, 'prefix') || ...
   ~ischar(opts.prefix) || ...
    numel(opts.prefix) > 1
    opts.prefix = '';
end
ppot = class(ppo);
if ischar(ppo)
    ppo = cellstr(ppo);
end
oppo = ppo;
ppoc = [];

% loop for multiple objects
if iscell(ppo)
    ppoc = ppo(:);
    
    % single object after all
    if numel(ppo) == 1
        
        % load object
        if ischar(ppo{1})
        end
        
        % was char
        if strcmpi(ppot, 'char')
            oppo = oppo{1};
        end
        
        % return
        return;
    end
    
    % otherwise
    try
        ppo = xff(ppoc);
    catch ne_eo;
        rethrow(ne_eo);
    end
end

% single (xff) object
if numel(ppo) == 1
    
    % type
    ppotype = lower(ppo.Filetype);
    switch (ppotype)
        
        % DMR
        case {'dmr'}
            
            % convert data
            ppo
            
        % FMR
        case {'fmr'}
            
        % HDR/NII
        case {'hdr'}
            
        % HEAD (BRIK/AFNI)
        case {'head'}
            
        % VTC
        case {'vtc'}
    end
    % prefix?
    if ~isempty(opts.prefix)
        
        % first 
    end
end

% iterate over each object
for oc = 1:numel(ppo)
end
