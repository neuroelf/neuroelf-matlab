function analyze_setendian(hdrfiles, entype, renfrom, rento)
% analyze_setendian  - set endian type for Analyze files
%
% FORMAT:       analyze_setendian(hdrfiles, entype [, renfrom, rento])
%
% Input fields:
%
%       hdrfiles    cell array with HDR filenames
%       entype      requested endian type (either 'ieee-be' or 'ieee-le')
%       renfrom     if given, rename from "pattern"
%       rento       if given, rename to "pattern"
%
% No output fields.
%
% Note: if from pattern contains * or + uses regexprep instead of strrep
% Note: if to pattern contains % uses sprintf with file number argument

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-08 2011, 10:18 PM EST
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

% argument check
if nargin < 2 || ...
   ~iscell(hdrfiles) || ...
    isempty(hdrfiles) || ...
   ~ischar(hdrfiles{1}) || ...
    isempty(hdrfiles{1}) || ...
   ~ischar(entype) || ...
   ~any(strcmpi(entype(:)', {'b', 'be', 'ieee-be', 'ieee-le', 'l', 'le'}))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
hdrfiles = hdrfiles(:);
if any(lower(entype(:)) == 'b')
    entype = 'ieee-be';
else
    entype = 'ieee-le';
end
entype = lower(entype(:)');
fn = numel(hdrfiles);
for fc = fn:-1:1
    if ~ischar(hdrfiles{fc}) || ...
        numel(hdrfiles{fc}) < 5 || ...
       ~strcmpi(hdrfiles{fc}(end-3:end), '.hdr') || ...
        exist(hdrfiles{fc}(:)', 'file') ~= 2
        hdrfiles(fc) = [];
    else
        hdrfiles{fc} = hdrfiles{fc}(:)';
    end
end
if isempty(hdrfiles)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin > 3 && ...
    ischar(renfrom) && ...
   ~isempty(renfrom) && ...
    ischar(rento)
    doren = true;
    renfrom = renfrom(:)';
    rento = rento(:)';
    if any(renfrom == '+' | renfrom == '*')
        useregx = true;
    else
        useregx = false;
    end
    if any(rento == '%')
        usespr = true;
    else
        usespr = false;
    end
else
    doren = false;
end

% iterate over files
for fc = 1:numel(hdrfiles)

    % try loading HDR file
    try
        hdr = [];
        hdr = xff(hdrfiles{fc});
        hdr.LoadVoxelData;
    catch ne_eo;
        if isxff(hdr)
            hdr.ClearObject;
        end
        error( ...
            'neuroelf:BadFileContent', ...
            'Bad HDR file content in %s (%s).', ...
            hdrfiles{fc}, ne_eo.message ...
        );
    end

    % continue if not rename and endian type the same
    if ~doren && ...
        strcmpi(entype, hdr.Endian)
        hdr.ClearObject;
        continue;
    end

    % set endian type
    hdr.Endian = entype;

    % perform rename if required
    if ~doren
        nname = hdr.FilenameOnDisk;
    else
        [oname{1:3}] = fileparts(hdr.FilenameOnDisk);
        if isempty(oname{1})
            oname{1} = '.';
        end
        if useregx
            if usespr
                nname = regexprep(oname{2}, renfrom, sprintf(rento, fc));
            else
                nname = regexprep(oname{2}, renfrom, rento);
            end
        else
            if usespr
                nname = strrep(oname{2}, renfrom, sprintf(rento, fc));
            else
                nname = strrep(oname{2}, renfrom, rento);
            end
        end
        nname = [oname{1} '/' nname oname{3}];
    end

    % save as
    hdr.SaveAs(nname);
    hdr.ClearObject;
end
