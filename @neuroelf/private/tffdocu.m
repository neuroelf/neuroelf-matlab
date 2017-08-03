function [varargout] = tffdocu
% TFF documentation (file format, usage)
% -------------------------------------------------------------------
%
% TFF (text file format) is a text based specification given to
% read text files (primarily with fixed content fields).
%
% The specification file itself must meet the following requirements
% to be recognized (and successfully parsed by tffparse()):
%
%  - empty lines are removed from the file
%  - all lines will be cut after any occurance of a hash mark (#)
%  - the first line of the file MUST contain the token 'TextFileFormat'
%  - the following tokens are recognized as single line options:
%    * CustomDelimiters:{[cc, ...], ...}  (character code sequences)
%    * FieldDelimiters:{[cc, ...], ...}   ('')
%    * LineDelimiters:{[cc, ...], ...}    ('')
%    * ArrayFormat:<sprintf format string for array element>
%    * BinaryIO:0|1 (do not read the file as ASCII/split cells)
%    * Extensions:<comma separated list of supported file extensions>
%    * FieldDelimCollapse:0|1 (collapse multiple field delimiters to one)
%    * FilenameMatch:<comma separated list of regexpi tokens>
%    * Filetype:<string used for file open dialogs and texts>
%    * ParagraphArrays:0|1 (automatically insert a blank line after array)
%    * SkipEmptyLines:0|1 (skip empty lines during read)
%  - the following tokens are used as section delimiters
%    * AfterReadCode:                     [opens a section until]
%    * EndAfterReadCode
%    * BeforeWriteCode:                   [opens a section until]
%    * EndBeforeWriteCode
%    * ListOfFields:<field delimiter>     [opens a section until]
%    * EndListOfFields
%    * Magic:<field delimiter>            [opens a section until
%    * EndMagic:
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
% text file contents. Conditional reads and writes can be made (to
% support different sub-formats or versions of the same file type),
% loops can be used and valid expressions calculated. The section
% itself is a formatted table which in its first row must contain at
% least the following header tokens (order irrelevant):
%
%  - type:     either of
%              - 'ARRAY': unnamed rectangular matrix with multiple values
%              - 'BLOOP': begin loop (name in varname, length in dim)
%              - 'ELOOP': end loop (name must match)
%              - 'EXPRE': evaluate an expression
%              - 'FIELD': read a mandatory field from file
%              - 'FLIST': read an optional field from file (out of a list)
%              - 'SKIPN': skip the next N rules (given in dim)
%              - 'WRTLN': write specific line to output file (in varname)
%              - 'XLOOP': exit loop (on condition, name must match)
%              lines with an empty type field are discarded
%  - cond:     conditional statement when to go/exit into a loop or
%              read/write a field, e.g. '@FileVer > 1 & @FileVer < 16'
%  - field:    name of field, as should appear in text file;
%              if left empty, copied from varname
%  - datatype: storage datatype (uint32, double, colorcode, ...)
%              custom datatypes, for which matching functions MUST exist,
%              can be given (string2DATATYPE, DATATYPE2string, and
%              array2DATATYPE, and DATATYPE2array), e.g. string2colorcode
%  - format:   sprintf output format of either field or array element
%  - dim:      dimension (1-D for BLOOP, FIELD, and FLIST; 2-D for ARRAY)
%  - default:  default value (in datatype syntax), if given, FIELD
%              content may be missing on writing calls to tffio(...)
%  - varname:  depending on type:
%              - 'BLOOP': loop variable name (used by $LOOPVAR)
%              - 'EXPRE': evaluated expression
%              - 'ARRAY', 'FIELD', 'FLIST': source/target variable
%
% Next to the built-in datatypes (fread, fwrite), one additional
% datatype is supported by tffio: string. It reads/writes until the end
% of the line. Multiple strings can be used in an Nx1 ARRAY.
%
% EXPREssions can re-use FIELD varnames (FIELD itself is able to use
% inline expressions). Variables are addressed via $VARNAME and,
% internally, resolved via the struct namevars.(...).
%
% FIELDs are read into both namevars.(...) and tffcont.(...). To set
% or alter content in the tffcont struct, the syntax @VARNAME can be
% used in EXPREssions.
%
% The following variables can be used as predefined variables:
%
% $TFFVERSION: version string of tffio, also copied to tffcont struct
% $TFFREAD:    if true, indicates that the file is being read from
% $TFFWRITE:   if true, indicates that the file is being written to
% $EXTENSION:  file extension of the file currently being read/written
% $FILENAME:   name of the file currently being read/written
%
% # Here is an example for a ListOfFields section:
% ListOfFields:!
% type !cond        !field       !datatype!format !dim    !default !varname
%
% # read file version
% FIELD!            !            !double  !%d     !1      !1       !FileVersion
%
% # get dimensions for writing
% EXPRE!$TFFWRITE   !!!!!!@DimX = size(@BinData, 1);
% EXPRE!$TFFWRITE   !!!!!!@DimY = size(@BinData, 2);
%
% # read/write dimensions
% FIELD!            !            !double  !%.0g   !1      !        !DimX
% FIELD!            !            !double  !%.0g   !1      !        !DimY
%
% # read/write data
% ARRAY!            !            !uint8   !%d     !@DimX, @DimY !  !BinData
% EndListOfFields
%
% For more extensive examples, please have a look at the TFF spec
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
%              - regexp:  regular expression, any 0x{CC} are replaced
%                         to allow custom characters, such as | and #
%              - regexpi: regular expression, case ignored
%
% # Here is an example for a Magic section
% Magic:|
% name         |range       |type    |magic
% PRT_name     |1, 512      |strfind |ProtocolName:
% EndMagic
%
% In practice this would mean that a file containing the string
% 'ProtocolName:' within the byte range [1, 512] would be recognized
% as, in this case, a BrainVoyager QX PRT (protocol) file.

% Version:  v0.9b
% Build:    10061507
% Date:     Jun-15 2010, 12:37 PM EST
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
help tffdocu;
