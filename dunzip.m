function M = dunzip(Z, zclass, zsize)
% DUNZIP - decompress DZIP output to recover original data
%
% USAGE:
% M = dzip(Z)
%
% VARIABLES:
% Z = compressed variable to decompress
% M = decompressed output
%
% NOTES: (1) The input variable Z is created by the DZIP function and
%            is a vector of type uint8
%        (2) The decompressed output will have the same data type and
%            dimensions as the original data provided to DZIP.
%        (3) See DZIP for other notes.
%        (4) Carefully tested, but no warranty; use at your own risk.
%        (5) Michael Kleder, Nov 2005
%
% see https://www.mathworks.com/matlabcentral/fileexchange/8899

classes = {'double','single','logical','char','int8','uint8', ...
    'int16','uint16','int32','uint32','int64','uint64'};

import com.mathworks.mlwidgets.io.InterruptibleStreamCopier
a=java.io.ByteArrayInputStream(Z);
b=java.util.zip.InflaterInputStream(a);
isc = InterruptibleStreamCopier.getInterruptibleStreamCopier;
c = java.io.ByteArrayOutputStream;
isc.copyStream(b,c);
Q=typecast(c.toByteArray,'uint8');
datastart = 1;
if nargin < 3 || ~ischar(zclass) || ~any(strcmpi(classes, zclass)) || ...
   ~isa(zsize, 'double') || numel(zsize) < 2 || size(zsize, 2) ~= numel(zsize) || ...
    any(isnan(zsize) | zsize < 1 | zsize ~= fix(zsize))
    if nargin == 3
        error('neuroelf:dunzip:badArgument', 'Bad zclass or zsize argument.');
    end
    cn = double(Q(1)); % class
    nd = double(Q(2)); % # dims
    zsize = typecast(Q(3:8*nd+2),'double')'; % size
    datastart = datastart + 2 + 8 * nd;
else
    cn = find(strcmpi(classes, zclass));
end
if datastart > 1
    Q=Q(datastart:end);
end
if cn == 3
    M  = logical(Q);
elseif cn == 4
    M = char(Q);
else
    M = typecast(Q,classes{cn});
end
try
    if any(isinf(zsize))
        nisize = zsize(~isinf(zsize));
        rsize = numel(M) / prod(nisize);
        if rsize ~= fix(rsize)
            return
        end
        zsize(isinf(zsize)) = rsize;
    end
    M=reshape(M,zsize);
end
return
