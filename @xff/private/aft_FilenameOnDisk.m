function fname = aft_FilenameOnDisk(xo, source)
% AFT::FilenameOnDisk  - returns the filename property
%
% FORMAT:       filename = obj.FilenameOnDisk;
%
% No input fields
%
% Output fields:
%
%       filename    last used filename (on load or save)
%
% TYPES: ALL

% Version:  v1.0
% Build:    16012414
% Date:     Jan-24 2016, 2:00 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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

% only valid for single file
if numel(xo) ~= 1 || ~xffisobject(xo, true)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% or of source
if nargin > 1 && numel(source) == 1 && ((islogical(source) && source) || ...
    (isa(source, 'double') && ~isinf(source) && ~isnan(source) && source ~= 0)) && ...
    isfield(xo.H, 'GZIPfile') && ischar(xo.H.GZIPfile) && ~isempty(xo.H.GZIPfile)
    if isa(source, 'double') && source == 2 && isfield(xo.H, 'GZIPext') && ...
        ischar(xo.H.GZIPext) && ~isempty(xo.H.GZIPext)
        fname = [xo.H.GZIPfile(:)', xo.H.GZIPext(:)'];
    else
        fname = xo.H.GZIPfile(:)';
    end

% get file name
else
    fname = xo.F;
end
