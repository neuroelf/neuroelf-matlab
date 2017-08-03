function [xfff, xffconf] = xffformats(fdir, xffconf)
% xff::_formats  - reading xff formats from folder
%
% THIS FUNCTION IS AN INTERNAL FUNCTION OF THE CLASS
%
% @xff
%
% AND IT SHOULD NEVER BE CALLED FROM OUTSIDE THE CLASS CODE

% Version:  v1.1
% Build:    16021216
% Date:     Feb-12 2016, 4:33 PM EST
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

% use global config
global ne_methods;

% argument check
if nargin < 2 || ~ischar(fdir) || isempty(fdir) || exist(fdir(:)', 'dir') ~= 7
    error('neuroelf:xff:badArgument', ...
        'xffformats requires a non-empty 1xN char folder argument.');
end

% make valid folder name
fdir = fdir(:)';
if fdir(end) ~= filesep
    fdir(end + 1) = filesep;
end

% read format files
try
    bff = dir([fdir '*.bff']);
    tff = dir([fdir '*.tff']);

    % initialize global counter, extension list, and magic tokens array
    bnf = 0;
    tnf = 0;
    bffspec = struct;
    tffspec = struct;
    ex = struct;
    gx = [];
catch xfferror
    error('neuroelf:xff:errorReadingDirContents', ...
        'Could not get contents of folder ''%s'' (%s).', fdir, xfferror.message);
end

% init output argument
xfff = struct;

% parsing found BFF files
for fc = 1:numel(bff)

    % reading singular format
    try
        drawnow;
        sf = ne_methods.bffparse([fdir bff(fc).name]);

        % store parsed info
        bnf = bnf + 1;
        if bnf == 1
            bffspec = sf;
            bffspec(numel(bff)).FFTYPE = 'BFF';
        else
            bffspec(bnf) = sf;
        end

        % add valid extensions -> singular format, and type config
        for ec = 1:length(sf.Extensions)
            lex = lower(sf.Extensions{ec});
            ex.(lex) = {bff(fc).name, bnf};
            xffconf.type.(lex) = struct;
            xffconf.update.(lex) = true;
        end

        % put magic tokens into main array
        for mc = 1:length(sf.Magic)
            if isempty(gx)
                gx = sf.Magic(mc);
            else
                gx(end+1) = sf.Magic(mc);
            end
        end
    catch xfferror
        warning('neuroelf:xff:badFileContents', ...
            'Bad file contents for ''%s'' (%s), please re-install.', ...
            bff(fc).name, xfferror.message);
    end
end

% parsing found TFF files
for fc = 1:numel(tff)

    % reading singular format
    try
        drawnow;
        sf = ne_methods.tffparse([fdir tff(fc).name]);

        % store parsed info
        tnf = tnf + 1;
        if tnf == 1
            tffspec = sf;
            tffspec(numel(tff)).FFTYPE = 'TFF';
        else
            tffspec(tnf) = sf;
        end

        % create empty settings struct for AFT
        xffconf.type.aft = struct;

        % add valid extensions -> singular format
        for ec = 1:length(sf.Extensions)
            lex = lower(sf.Extensions{ec});
            ex.(lex) = {tff(fc).name, tnf};
            xffconf.type.(lex) = struct;
            xffconf.update.(lex) = true;
        end

        % put magic tokens into main array
        for mc = 1:length(sf.Magic)
            if isempty(gx)
                gx = sf.Magic(mc);
            else
                gx(end+1) = sf.Magic(mc);
            end
        end
    catch xfferror
        warning('neuroelf:xff:badFileContents', ...
            'Bad file contents for ''%s'' (%s), please re-install.', ...
            tff(fc).name, xfferror.message);
    end
end

% set output fields
xfff.bff = bffspec;
xfff.tff = tffspec;
xfff.extensions = ex;
xfff.magic = gx;
