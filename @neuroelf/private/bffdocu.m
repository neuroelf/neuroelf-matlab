function [varargout] = bffdocu
% BFF documentation (file format, usage)
% -------------------------------------------------------------------
%
% BFF (binary file format) is a text based specification given to
% read binary files (primarily with fixed sized content).
%
% The specification file itself must meet the following requirements
% to be recognized (and successfully parsed by bffparse()):
%
%  - empty lines are removed from the file
%  - all lines will be cut after any occurance of a hash mark (#)
%  - the first line of the file MUST contain the token 'BinaryFileFormat'
%  - the following tokens are recognized as single line options:
%    * EncodingSyntax:<either 'native', 'ieee-le', or 'ieee-be'>
%    * Extensions:<comma separated list of supported file extensions>
%    * FilenameMatch:<comma separated list of regexpi tokens>
%    * Filetype:<string used for file open dialogs and texts>
%    * TransIOSize:<string that must evaluate to a number of bytes>
%  - the following tokens are used as section delimiters
%    * AfterReadCode:                     [opens a section until]
%    * EndAfterReadCode
%    * BeforeWriteCode:                   [opens a section until]
%    * EndBeforeWriteCode
%    * ListOfFields:<field delimiter>     [opens a section until]
%    * EndListOfFields
%    * Magic:<field delimiter>            [opens a section until]
%    * EndMagic
%
% ----
%
% The AfterReadCode and BeforeWriteCode sections may contain code
% that is executed either after reading the file or before writing
% the file to disk (even before a file open is attempted in the
% latter case).
%
% Both sections may contain ($|@)([a-z][a-z_0-9]*)/i variable
% references, which will be resolved correctly. Make sure that the
% code does NOT contain hash marks (other than to comment out things
% that is!)
%
% ----
%
% The ListOfFields section describes the rules to read and write the
% binary file contents. Conditional reads and writes can be made (to
% support different sub-formats or versions of the same file type),
% loops can be used and valid expressions calculated. The section
% itself is a formatted table which in its first row must contain at
% least the following header tokens (order irrelevant):
%
%  - type:     either of
%              - 'BLOOP': begin loop (name in varname, length in dim)
%              - 'ELOOP': end loop (name must match)
%              - 'EXPRE': evaluate an expression
%              - 'FIELD': read a field from file
%              - 'SKIPN': skip the next N rules (given in dim)
%              - 'XLOOP': exit loop (on condition, name must match)
%              lines with an empty type field are discarded
%  - cond:     conditional statement when to go/exit into a loop or
%              read/write a field, e.g. '@FileVer > 1 & @FileVer < 16'
%  - disktype: reading datatype (char, uint16, int32, single, ...)
%  - datatype: storage datatype (uint32, double, colorcode, ...)
%              custom datatypes can be given, for which two functions
%              MUST exist (DISKTYPE2DATATYPE and DATATYPE2DISKTYPE,
%              e.g. uint322colorcode and colorcode2uint32)
%  - dim:      dimension (1-D for BLOOP, N-D for FIELD)
%  - default:  default value (in datatype syntax), if given, FIELD
%              content may be missing on writing calls to bffio(...)
%  - varname:  depending on type:
%              - 'BLOOP': loop variable name (used by $LOOPVAR)
%              - 'EXPRE': evaluated expression
%              - 'FIELD': source/target variable (postfixed to struct)
%
% Next to the built-in datatypes (fread, fwrite), one additional
% datatype is supported by bffio: cstring. It reads until the next
% 0x00 value is encountered in the stream, and upon writing adds this
% end-of-string signature into the file. Multiple strings can be used
% by giving a non singleton dimension.
%
% EXPREssions can re-use FIELD varnames (FIELD itself is able to use
% inline expressions). Variables are addressed via $VARNAME and,
% internally, resolved via the struct namevars.(...).
%
% FIELDs are read into both namevars.(...) and bffcont.(...). To set
% or alter content in the bffcont struct, the syntax @VARNAME can be
% used in EXPREssions.
%
% The following variables can be used as predefined variables:
%
% $BFFVERSION: version string of bffio, also copied to bffcont struct
% $BFFREAD:    if true, indicates that the file is being read from
% $BFFWRITE:   if true, indicates that the file is being written to
% $EXTENSION:  file extension of the file currently being read/written
% $FILENAME:   name of the file currently being read/written
% $FILESIZE:   size of the file (only available in read mode)
% $HEADERONLY: if true, indicates that only header fields are to be read
%
% # Here is an example for a ListOfFields section:
% ListOfFields:!
% type !cond        !disktype!datatype!dim       !default !varname
%
% # read file version
% FIELD!            !uint16  !double  !1, 1      !1       !FVersion
%
% # get dimensions for writing
% EXPRE!$BFFWRITE   !!!!!@DimX = size(@BinData, 1);
% EXPRE!$BFFWRITE   !!!!!@DimY = size(@BinData, 2);
%
% # read/write dimensions
% FIELD!            !uint32  !double  !1, 1      !        !DimX
% FIELD!            !uint32  !double  !1, 1      !        !DimY
%
% # read/write data
% FIELD!            !uint8   !uint8   !@DimX, @DimY !     !BinData
% EndListOfFields
%
% For more extensive examples, please have a look at the BFF spec
% files coming with this toolbox.
%
% ----
%
% The Magic section describes, if possible for the format, how it can
% be detected from its content. The section itself is a formatted table
% which in its first row must contain at least the following header
% tokens (order irrelevant):
%  - name:     unique token to go with this format
%  - range:    1x2 double (integer) array to look for in the file
%              negative numbers are considered as reverse from EOF
%  - type:     either of 'strfind', 'hex', 'regexp'
%  - magic:    pattern to match to accept content;
%              - strfind: simple string
%              - hex:     packed content, syntax: 00,08,00,10[,...]
%              - hex16:   packed content, syntax: 0800,0010[,...]
%              - hex32:   packed content, syntax: 08000010[,...]
%              - regexp:  regular expression, any 0x{CC} are replaced
%                         to allow custom characters, such as | and #
%              - regexpi: regular expression, case ignored
%              - strfind: do a simple strfind with range and pattern
%
% # Here is an example for a Magic section
% Magic:|
% name         |range       |type    |magic
% BMP_bmType   |1, 2        |strfind |BM
% EndMagic
%
% In practice this would mean that a file containing the string 'BM'
% within the byte range [1, 2] would be recognized as, in this case, a
% Windows BITMAP (BMP) file.

% Version:  v0.9a
% Build:    10051716
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

varargout = cell(1, nargout);
help bffdocu;
