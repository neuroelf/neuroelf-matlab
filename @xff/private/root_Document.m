function xo = root_Document(xo, dspec, fload)
% ROOT::Document  - get one "Document" (VB-Style interface)
%
% FORMAT:       object = xff.Document(dspec);
%
% Input fields:
%
%       dspec       either numbered object or (partial) filename
%
% Output fields:
%
%       object      found object (otherwise: error)
%
% Using: findfirst.

% Version:  v1.1
% Build:    16053018
% Date:     May-30 2016, 6:07 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% neuroelf library and factory
global ne_methods xffsngl;
findfirst = ne_methods.findfirst;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'root') || ...
   (~isa(dspec, 'double') && ~ischar(dspec)) || isempty(dspec)
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if nargin < 3 || ~islogical(fload) || numel(fload) ~= 1
    fload = false;
end

% get available objects (without ROOT!)
o = xff(0, 'objects');

% no objects loaded yet
if isempty(o)

    % try to load anyway
    if ischar(dspec) && ~isempty(dspec) && exist(dspec, 'file') == 2 && fload
        try
            xo = xff(dspec(:)');
            return;
        catch xfferror;
            neuroelf_lasterr(xfferror);
        end
    end
    error('neuroelf:xff:documentNotFount', 'Document not found.');
end

% for strings, look for xffID, filename and partial match
if ischar(dspec)
    dspec = dspec(:)';
    fm = [];
    if numel(dspec) == 24 && all(dspec >= '0' & upper(dspec) <= 'F')
        ons = xffsngl.OBJS(2:end, 3);
        fm = findfirst(strcmpi(dspec, ons));
    end
    if isempty(fm)
        ons = {o(:).F};
        fm = findfirst(strcmpi(dspec, ons));
    end
    if isempty(fm)
        for c = 1:numel(ons)
            [p, f, e] = fileparts(ons{c});
            ons{c} = [f, e];
        end
        [p, f, e] = fileparts(dspec);
        f = [f, e];
        fm = findfirst(strcmpi(f, ons));
    end
    if ~isempty(fm)
        xo = o(fm);
    elseif ~isempty(dspec) && exist(dspec, 'file') == 2 && fload
        try
            xo = xff(dspec(:)');
            return;
        catch xfferror;
            neuroelf_lasterr(xfferror);
        end
    else
        error('neuroelf:xff:documentNotFount', 'Document not found.');
    end

% invalid double index lookup
elseif numel(dspec) ~= 1 || isinf(dspec) || isnan(dspec) || dspec < 1 || dspec > numel(o)
    xo = o([]);

% valid index
else

    % make object
    xo = o(round(dspec));
end

% error out if empty
if isempty(xo)
    error('neuroelf:xff:documentNotFount', 'Document not found.');
end
