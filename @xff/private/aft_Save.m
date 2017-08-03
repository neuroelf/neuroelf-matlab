function xo = aft_Save(xo)
% AnyFileType::Save  - saves any xff object back to disk
%
% FORMAT:       object.Save;
%
% No input / output fields.
%
% TYPES: ALL
%
% Using: bffio, tffio.

% Version:  v1.1
% Build:    16040820
% Date:     Apr-08 2016, 8:01 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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

% global neuroelf methods
global ne_methods;

% only valid for single file
if numel(xo) ~= 1 || ~xffisobject(xo, true) || xo.L(1) == 'X' || ~isfield(xo.S, 'FFTYPE')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% check filename
if isempty(xo.F)
    error('neuroelf:xff:badFilename', 'File not yet saved. Use SaveAs method instead.');
end

% don't allow volume marker
if ~isempty(regexpi(xo.F, ',\d+$'))
    error('neuroelf:xff:badArgument', 'Saving of sub-volumes not permitted.');
end

% what to do
try
    switch (lower(xo.S.FFTYPE))
        case 'bff'
            xo.C = ne_methods.bffio(xo.F, xo.S, xo.C);
        case 'tff'
            [xo.C, xo.F] = ne_methods.tffio(xo.F, xo.S, xo.C);
        otherwise
            error('neuroelf:xff:invalidFileType', 'Type not recognized (?FF): %s.', xo.S.FFTYPE);
    end
catch xfferror
    error('neuroelf:xff:errorSavingFile', 'Error saving file %s: %s.', xo.F, xfferror.message);
end

% pack
if isfield(xo.H, 'GZIPext') && ischar(xo.H.GZIPext) && strcmpi(xo.H.GZIPext, '.gz') && ...
    isfield(xo.H, 'GZIPfile') && ischar(xo.H.GZIPfile) && ~isempty(xo.H.GZIPfile)
    try
        gzip(xo.F);
        [cps, cpm, cpi] = copyfile([xo.F '.gz'], [xo.H.GZIPfile xo.H.GZIPext]);
        if cps ~= 1
            error(cpi, cpm);
        end
    catch xfferror
        rethrow(xfferror);
    end
end

% then see if RunTimeVars are to be saved as well
if isfield(xo.C.RunTimeVars, 'AutoSave') && islogical(xo.C.RunTimeVars.AutoSave) && ...
    numel(xo.C.RunTimeVars.AutoSave) == 1 && xo.C.RunTimeVars.AutoSave

    % try automatic saving
    try
        aft_SaveRunTimeVars(xo);
    catch xfferror
        warning('neuroelf:xff:errorSavingFile', 'Error saving RunTimeVars file: %s.', xfferror.message);
    end
end
