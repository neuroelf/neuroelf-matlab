function [xo, sorder] = fmr_ReorderSlices(xo, sorder, sgap)
% FMR::ReorderSlices  - reorders the slices in an FMR
%
% FORMAT:       [fmr, bwsorder] = fmr.ReorderSlices(sorder [, sgap, snum])
%
% Input fields:
%
%       sorder      either an 1xS array with new order
%                   or 1x1 value: 1 for forward, -1 for backwards, needs
%       sgap        slice enumeration gap, default: round(sqrt(NrOfSlices))
%
% Output fields:
%
%       fmr         fmr with altered slice order
%       bwsorder    backwards sorting order
%
% Note: this function can be used when none of the default slice orderings
%       for slice scan-time correction in BrainVoyager QX matches the true
%       order of the acquired scan.
%
%       for DataStorageFormat 1, this function simply renames the slice
%       files (and reloads the FMR slices); for DataStorageFormat 2, the
%       slices are reordered within the STC file!
%
% Using: renamefile.

% Version:  v1.1
% Build:    16020311
% Date:     Feb-03 2016, 11:27 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% neuroelf library
global ne_methods;
renamefile = ne_methods.renamefile;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'fmr') || ...
   ~isa(sorder, 'double') || (numel(sorder) ~= 1 && numel(sorder) ~= max(size(sorder))) || ...
    any(isinf(sorder) | isnan(sorder) | sorder < -1 | sorder ~= fix(sorder)) || ...
    numel(sorder) ~= numel(unique(sorder))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
sorder = sorder(:)';

% get FMR content
bc = xo.C;

% check sorder field for details
if numel(sorder) ~= 1 && numel(sorder) ~= bc.NrOfSlices
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if numel(sorder) == 1
    if nargin < 3 || ~isa(sgap, 'double') || numel(sgap) ~= 1 || isnan(sgap) || ...
        sgap < 1 || sgap >= bc.NrOfSlices || sgap ~= fix(sgap)
        sgap = round(sqrt(bc.NrOfSlices));
    end
    if sorder < 0
        backw = true;
    else
        backw = false;
    end
    sorder = 1:sgap:bc.NrOfSlices;
    for gc = 2:sgap
        sorder = [sorder, gc:sgap:bc.NrOfSlices];
    end
    if backw
        [orderi{1:2}] = sort(sorder);
        sorder = orderi{2};
    end
end

% for slices in separate files
if bc.DataStorageFormat == 1

    % change to FMR folder
    opwd = pwd;
    npwd = fileparts(xo.F);
    try
        cd(npwd);
    catch xfferror
        neuroelf_lasterr(xfferror);
        error('neuroelf:xff:cdError', 'Error changing directory to %s.', npwd);
    end

    % check whether slice files exist
    for sc = 1:bc.NrOfSlices
        if exist(sprintf('%s-%d.stc', bc.Prefix, sc), 'file') ~= 2
            try
                cd(opwd)
            catch xfferror
                neuroelf_lasterr(xfferror);
            end
            error('neuroelf:xff:referencedFileNotFound', ...
                'Referenced STC file of slice %d not found.', sc);
        end
    end

    % rename slice files (two passes!)
    for sc = 1:bc.NrOfSlices
        try
            renamefile(sprintf('./%s-%d.stc', bc.Prefix, sorder(sc)), ...
                sprintf('./R_%s-%d.stc', bc.Prefix, sc));
        catch xfferror
            error('neuroelf:xff:renameError', ...
                'Error renaming slice %d to %d: %s.', sorder(sc), sc, xfferror.message);
        end
    end
    for sc = 1:bc.NrOfSlices
        try
            renamefile(sprintf('./R_%s-%d.stc', bc.Prefix, sc), ...
                sprintf('./%s-%d.stc', bc.Prefix, sc));
        catch xfferror
            error('neuroelf:xff:renameError', ...
                'Error final-renaming slice %d: %s.', sorder(sc), sc, xfferror.message);
        end
    end

    % reload STC (transio) from data
    fmr_LoadSTC(xo);
    bc = xo.C;

    % change back to original folder
    try
        cd(opwd);
    catch xfferror
        neuroelf_lasterr(xfferror);
    end

% for single-STC file
else
    if istransio(bc.Slice.STCData)
        STCData = resolve(bc.Slice.STCData);
        STCData = STCData(:, :, :, sorder);
        bc.Slice.STCData(:, :, :, :) = STCData;
    else
        bc.Slice.STCData = bc.Slice.STCData(:, :, :, sorder);
        stcf = [fileparts(xo.F) '/' bc.Prefix '.stc'];
        try
            stc = 0;
            stc = fopen(stcf, 'w', 'ieee-le');
            if stc < 1
                error('ERROR_WRITING_STC');
            end
            if bc.DataType == 1
                fwrite(stc, bc.Slice.STCData(:), 'uint16');
            else
                fwrite(stc, bc.Slice.STCData(:), 'single');
            end
            fclose(stc);
        catch xfferror
            if stc > 0
                fclose(stc);
            end
            error('neuroelf:xff:ErrorWritingFile', ...
                'Error writing STC file %s.stc: %s', bc.Prefix, xfferror.message);
        end
        xo.C = bc;
    end
end

% reorder AMR as well?
amrfile = bc.LoadAMRFile;
if ischar(amrfile) && ~isempty(amrfile) && exist(amrfile, 'file') == 2
    try
        amr = [];
        amr = xff(amrfile);
        abc = amr.C;
        abc.Slice = abc.Slice(sorder);
        amr.C = abc;
        aft_Save(amr);
    catch xfferror
        delete(amr);
        warning('neuroelf:xff:internalError', ...
            'Error applying reordering to linked AMR: %s.', xfferror.message);
    end
end

% reverse sorting as output argument
[orderi{1:2}] = sort(sorder);
sorder = orderi{2};
