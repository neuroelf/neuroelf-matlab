function readhrv(filename, targetdir, onsets, seglen, options)
% readhrv  - reads a heart rate variability text file
%
% FORMAT:       readhrv(filename, targetdir, onsets, seglen [, options])
%
% Input fields:
%
%       filename    filename of HRV text file
%       targetdir   path to output files
%       onsets      1xN double array in seconds
%       seglen      1x1 or 1xN double with segment length (seconds)
%       options     1x1 struct with optional fields
%        .columns   1xN double array, column selection
%                   (default: 1..N, all columns in source order)
%        .decimals  1x1 double for number of decimals in output
%                   (default: 8, only used on .columns or .factor)
%	 .factor    1xN double array to multiply columns
%                   (default: 1 for every column detected)
%        .freq      frequency of data (default: 1000Hz)
%        .opat      output file pattern (default: seg%02d.txt)
%        .skip      1x1 double, skip N lines (default: 0)
%        .srcld     source line delimiter (default: auto, means
%                   detection of either LfCr, CrLf, Lf, Cr)
%        .trgld     target line delimiter (default: CrLf)
%        .writesn   [false|true], write scientific notation
%                   (default: true, required for .columns
%                   and/or .factor)
%
% See also: asciiread, asciiwrite.

% Version:  v0.9a
% Build:    10051716
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

% argument check
if nargin < 4 || ...
   ~ischar(filename) || ...
   ~ischar(targetdir) || ...
   ~isa(onsets, 'double') || ...
   ~isa(seglen, 'double') || ...
    exist(filename, 'file') ~= 2 || ...
    exist(targetdir, 'dir') ~= 7 || ...
    isempty(onsets) || ...
    isempty(seglen)
    error( ...
        'neuroelf:Badarguments', ...
        'readhrv requires four arguments. See ''help %s''.', ...
        mfilename ...
    );
end

if targetdir(end) ~= filesep
    targetdir = [targetdir filesep];
end

onsets = onsets(:)';
seglen = seglen(:)';
if length(seglen) == 1
    seglen = seglen * ones(size(onsets));
elseif length(seglen) < length(onsets)
    seglen = seglen(1) * ones(size(onsets));
else
    seglen = seglen(1:length(onsets));
end

defopt = { ...
    'columns',	'double',	'nonempty',	-1; ...
    'decimals',	'double',	'nonempty',	8; ...
    'factor',	'double',	'nonempty',	1; ...
    'freq',	'double',	'nonempty',	1000; ...
    'opat',	'char',		'nonempty',	'seg%02d.txt'; ...
    'skip',	'double',	'nonempty',     0; ...
    'srcld',	'char',		'',		''; ...
    'trgld',	'char',		'nonempty'	char([13, 10]); ...
    'writesn',	'logical',	'nonempty',	true ...
};

if nargin > 4 && ...
    isstruct(options) && ...
   ~isempty(options)
    options = checkstruct(options(1), defopt);
else
    options = checkstruct(struct, defopt);
end

options.columns = fix(options.columns(:)');
if any(options.columns<1)
    selcols = false;
else
    selcols = true;
end
options.decimals = fix(options.decimals(1));
if options.decimals < 0 || ...
    options.decimals > 12 || ...
    isinf(options.decimals) || ...
    isnan(options.decimals)
    options.decimals = 8;
end
options.factor = options.factor(:)';
if any(options.factor~=1)
    usefact = true;
else
    usefact = false;
end
options.freq = fix(options.freq(1));
if options.freq < 1
    options.freq = 1;
end
if options.freq > 10000
    options.freq = 10000;
end
if isinf(options.freq) || ...
    isnan(options.freq)
    options.freq = 1000;
end

if ~any(options.opat == '%')
    options.opat = 'seg%02d.txt';
end

options.skip = fix(options.skip(1));
if options.skip < 0 || ...
    isnan(options.skip) || ...
    isinf(options.skip)
    options.skip = 0;
end

disp('readhrv - heart rate variability file splitter');
disp(' ');
disp([' -> reading input file: ' filename]);
ar = asciiread(filename);

if isempty(options.srcld)
    arp = ar(1:4096);
    LfCr = strfind(arp, char([10, 13]));
    CrLf = strfind(arp, char([13, 10]));
    Lf = strfind(arp, char(10));
    Cr = strfind(arp, char(13));
    if ~isempty(LfCr)
        if isempty(CrLf)
            options.srcld = char([10, 13]);
        else
            if length(LfCr) > length(CrLf)
                options.srcld = char([10, 13]);
            else
                options.srcld = char([13, 10]);
            end
        end
    elseif ~isempty(CrLf)
        options.srcld = char([13, 10]);
    elseif ~isempty(Lf)
        options.srcld = char(10);
    elseif ~isempty(Cr)
        options.srcld = char(13);
    else
        error( ...
            'neuroelf:AutoDetectFailed', ...
            'Could not determine line delimiter in source.' ...
        );
    end
end

options.writesn = options.writesn(1);
if options.writesn && ...
   ~any(arp == ',')
    options.writesn = false;
end

sll = length(options.srcld);

disp([' -> finding line delimiters: ' any2ascii(double(options.srcld))]);
if sll == 1
    ap = [1 (find(ar==options.srcld)+1)];
else
    ap = [1 (strfind(ar, options.srcld)+sll)];
end
al = length(ap);

onsets = fix(onsets * options.freq) + 1 + options.skip;
seglen = fix(seglen * options.freq);
ofsets = onsets + seglen;

for c = length(ofsets):-1:1
    if ofsets(c) > al
        warning( ...
            'neuroelf:BadArgument', ...
            'Bad onset/seglen combination: #%d, [%d ... %d]', ...
            c, onsets(c), ofsets(c) ...
        );
        onsets(c) = [];
        seglen(c) = [];
        onsets(c) = [];
    end
end

if isempty(onsets)
    error( ...
        'neuroelf:BadArgument', ...
        'No valid segments requested.' ...
    );
end

disp(' -> splitting into segments...')
for c = 1:length(ofsets)
    as = ar(ap(onsets(c)):(ap(ofsets(c))-(sll+1)));
    if options.writesn
        as = strrep(as, ',', '.');
    end
    tfilen = [targetdir sprintf(options.opat, c)];
    asciiwrite(tfilen, strrep(as, options.srcld, options.trgld));
    disp(sprintf('    - segment %02d/%02d', c, length(ofsets)));
    if options.writesn && ...
       (selcols || usefact)
        ocols  = options.columns;
        ofact  = options.factor;
        oftld  = options.trgld;
        tfilec = load(tfilen);
        tfiles = size(tfilec);
        tcols  = tfiles(2);
        if ~any(ocols<1 | ocols>tcols)
            tfilenc = zeros([tfiles(1),length(ocols)]);
            for sc = 1:length(ocols)
                tfilenc(:,sc) = tfilec(:,ocols(sc));
            end
            tfilec = tfilenc;
            clear tfilenc;
        end
        if usefact && ...
            length(ofact) == size(tfilec, 2)
            for fc = 1:length(ofact)
                tfilec(:,fc) = tfilec(:,fc) .* ofact(fc);
            end
        end
        ofrm = {['%.' num2str(options.decimals) 'f']};
        ofrm = sprintf('%s\\t',ofrm{ones(1,size(tfilec,2))});
        ofrm(end) = 'n';
        fp = fopen(tfilen, 'wb');
        if (fp>0)
            numrows = size(tfilec,1);
            for rc = 1:1000:numrows
                fwrstr = sprintf(ofrm, tfilec(rc:min(rc+999,numrows),:)');
                if ~strcmp(oftld, sprintf('\n'))
                    fwrstr = strrep(fwrstr, sprintf('\n'), oftld);
                end
                fwrite(fp, fwrstr, 'char');
            end
            fclose(fp);
        end
    end
end
