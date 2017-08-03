function otype = ostype(tag)
% ostype  - returns information about the currently used OS
%
% FORMAT:       otype = ostype([tag])
%
% Input fields:
%
%       tag         optionally requesting a subfield of the struct
%
% Output fields:
%
%       otype       1x1 struct with fields
%        .computer  copy of computer()
%        .filesep   copy of filesep()
%        .ispc      copy of ispc()
%        .machine   either 'MAC', 'UNIX', or 'WIN'
%
% See also COMPUTER, FILESEP, ISPC

% Version:  v0.9a
% Build:    11043013
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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

% declare and (if needed) init persitently stored variable
persistent osinfo;
if isempty(osinfo) || ...
   ~isstruct(osinfo) || ...
   ~osinfo.is_initialized

    % initialize persistent variable
    osinfo = struct;
    osinfo.is_initialized = false;

    % get values
    osinfo.computer = upper(computer);
    osinfo.filesep = filesep;
    osinfo.ispc = ispc;

    if ~isempty(strfind(osinfo.computer, 'WIN')) && ...
        numel(strfind(osinfo.computer, 'WIN')) == 1 && ...
        strfind(osinfo.computer, 'WIN') <= 3
        osinfo.machine = 'WIN';
    elseif ~isempty(strfind(osinfo.computer, 'MAC'))
        osinfo.machine = 'MAC';
    else
        osinfo.machine = 'UNIX';
    end

    % sanity checks
    if (strcmp(osinfo.machine, 'UNIX') && filesep ~= '/') || ...
       (strcmp(osinfo.machine, 'WIN')  && filesep ~= '\')
        error( ...
            'neuroelf:InternalError', ...
            'Initialization failed. Bad filesep for machine type %s', ...
            osinfo.machine ...
        );
    end

    % set initialized to true
    osinfo.is_initialized = true;
end

% argument check
if nargin < 1 || ...
   ~ischar(tag) || ...
    isempty(tag)
    tag = 'full';
else
    tag = lower(tag(:));
end

switch (tag)
    case {'filesep'}
        otype = osinfo.filesep;
    case {'machine'}
        otype = osinfo.machine;
    otherwise
        otype = osinfo;
end
