function wavwrite2(y,Fs,nbits,wavefile)
%WAVWRITE Write Microsoft WAVE (".wav") sound file.
%   WAVWRITE(Y,WAVEFILE) writes a WAVE file specified by the
%   string WAVEFILE.  The data should be arranged with one channel
%   per column.  Amplitude values outside the range [-1,+1] are
%   clipped prior to writing.
%
%   WAVWRITE(Y,FS,WAVEFILE) specifies the sample rate FS of the data
%   in Hertz.
%
%   WAVWRITE(Y,FS,NBITS,WAVEFILE) forces an NBITS-bit file format to
%   be written, where NBITS<=16.
%
%   Supports multi-channel 8- or 16-bit WAVE data.
%
%   See also WAVREAD, AUWRITE.

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.6 $  $Date: 1997/11/21 23:24:12 $

%   D. Orofino, 11/95
%   revised by S.Frost 2/00 to change ftype/fcreator

% Get user default preferences:
Fs_pref = 8000; nbits_pref = 16;

% Parse inputs:
if nargin<2 | nargin>4,
  error('Incorrect number of input arguments.');
elseif nargin<3,
  wavefile=Fs;
  Fs=Fs_pref;
  nbits=nbits_pref;
elseif nargin<4,
  wavefile=nbits;
  nbits=nbits_pref;
end

% Open file for output:
if ~isstr(wavefile),
  error('wavefile must be a string.');
end
if isempty(findstr(wavefile,'.')),
  wavefile=[wavefile '.wav'];
end
fid=fopen(wavefile,'wb','l');  % Little-endian
if (fid==-1),
  error('Can''t open WAVE file for output.');
end

% If input is a vector, force it to be a column:
if ndims(y)>2,
  error('Data array must have 1- or 2-dimensions, only.');
end
[samples,channels]=size(y);
if samples==1,
   y=y(:); 
   [samples,channels]=size(y);
end

% Clip data to normalized range [-1,+1]:
i=find(abs(y)>1);
if ~isempty(i),
  y(i)=sign(y(i));
  warning('Data clipped during write to file.');
end

% # bytes per sample to write
bytes_per_sample = ceil(nbits/8);
total_samples = samples*channels;
total_bytes = total_samples * bytes_per_sample;

% Determine number of bytes in RIFF chunk
% (not including pad bytes, if needed):
% ----------------------------------
%  'RIFF'           4 bytes
%  size             4 bytes (ulong)
%  'WAVE'           4 bytes
%  'fmt '           4 bytes
%  size             4 bytes (ulong)
% <wave-format>     14 bytes
% <format_specific> 2 bytes (PCM)
%  'data'           4 bytes
%  size             4 bytes (ulong)
% <wave-data>       N bytes
% ----------------------------------
riff_cksize = 36+total_bytes;   % Don't include 'RIFF' or its size field
fmt_cksize = 16;                % Don't include 'fmt ' or its size field
data_cksize = total_bytes;      % Don't include 'data' or its size field

% Determine pad bytes:
data_pad = rem(data_cksize,2);
riff_cksize = riff_cksize + data_pad; % + fmt_pad, always 0

% Write RIFF chunk:
ck=[]; ck.fid=fid; ck.Size=riff_cksize; ck.ID='RIFF';
write_ckinfo(ck);

% Write WAVE:
ck.ID='WAVE';
write_ckinfo(ck,1);

% Write <fmt-ck>:
ck.ID='fmt ';
ck.Size=fmt_cksize;
write_ckinfo(ck);

% Write <wave-format>:
fmt.wFormatTag      = 1;            % Data encoding format = PCM
fmt.nChannels       = channels;     % Number of channels
fmt.nSamplesPerSec  = Fs;           % Samples per second
fmt.nAvgBytesPerSec = channels*bytes_per_sample*Fs; % Avg transfer rate
fmt.nBlockAlign     = channels*bytes_per_sample;    % Block alignment
fmt.nBitsPerSample  = nbits;        % standard <PCM-format-specific> info
status=write_wavefmt(fid,fmt);

% Write <data-ck>:
ck.ID='data';
ck.Size=data_cksize;
write_ckinfo(ck);

% Write <wave-data>, and its pad byte if needed:
status=write_wavedat(fid,fmt,y);
fclose(fid);    % Close file

% change ftype/fcreator
filetype(wavefile,'WAVE','Nqst');

% end of wavwrite()

% ------------------------------------------------------------------------
% Private functions:
% ------------------------------------------------------------------------

% WRITE_CKINFO: Writes next RIFF chunk, but not the chunk data.
%   If optional sflg is set to nonzero, write SUBchunk info instead.
%   Expects an open FID pointing to first byte of chunk header,
%   and a chunk structure.
%   ck.fid, ck.ID, ck.Size, ck.Data
function status=write_ckinfo(ck,sflg)
status=0;
if nargin<2, sflg=0; end
cnt = fwrite(ck.fid, ck.ID, 'char');
if cnt~=4, status=-1; end % Error condition
if ~sflg,
  % Write chunk size (skip if subchunk):
  cnt = fwrite(ck.fid, ck.Size, 'ulong');
  if cnt~=1, status=-1; end % Error condition
end
return;

% WRITE_WAVEFMT: Write WAVE format chunk.
%   Assumes fid points to the wave-format subchunk.
%   Requires chunk structure to be passed, indicating
%   the length of the chunk.
function status=write_wavefmt(fid,fmt)
status=0;

% Create <wave-format> data:
fwrite(fid, fmt.wFormatTag, 'ushort');
fwrite(fid, fmt.nChannels, 'ushort');
fwrite(fid, fmt.nSamplesPerSec, 'ulong');
fwrite(fid, fmt.nAvgBytesPerSec, 'ulong');
cnt=fwrite(fid, fmt.nBlockAlign, 'ushort');
if cnt~=1, status=-1; return; end   % Error condition

% Write format-specific info:
if fmt.wFormatTag==1,
  % Write standard <PCM-format-specific> info:
  cnt = fwrite(fid, fmt.nBitsPerSample, 'ushort');
  if cnt~=1, status=-1; end % Error condition
else
  error('Unknown data format.');
end

return;

% WRITE_WAVEDAT: Write WAVE data chunk
%   Assumes fid points to the wave-data chunk
%   Requires <wave-format> structure to be passed.
function status=write_wavedat(fid,fmt,data)
status=0;

if fmt.wFormatTag==1,
  % PCM Format:
  % Determine # bytes/sample - format requires rounding
  %  to next integer number of bytes:
  BytesPerSample = ceil(fmt.nBitsPerSample/8);
  if BytesPerSample==1,
    dtype='uchar'; % unsigned 8-bit
    % Scale data according to bits/samples: [-1,+1] -> [0,255]
    data=round((data+1) * 255/2);
  elseif BytesPerSample==2,
    dtype='short'; % signed 16-bit
    % Scale data according to bits/samples: [-1,+1] -> [-32768,+32767]
    data=round((data+1)*65535/2) - 32768;
  else
    error('Cannot write WAVE files with more than 16 bits/sample.');
  end

  % Write data, one row at a time (one sample from each channel):
  [samples,channels]=size(data);
  total_samples = samples*channels;
  cnt = fwrite(fid, reshape(data',total_samples,1), dtype);

  if cnt~=total_samples, status=-1; return; end % Error condition

  % Determine if a pad-byte is appended to data chunk:
  if rem(total_samples*BytesPerSample,2), fwrite(fid,0,'uchar'); end


else
  % Unknown wave-format for data.
  error('Unknown data format.');
end
return;

% end of wavwrite.m


