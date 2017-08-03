function mstr = mstring(xo, varargin)
% xfigure::mstring  - manipulate text of a listbox or dropdown
%
% FORMAT:       [subset] = mstring(hControl [,pos [,setstr [,insert]]])
%
% Input fields:
%
%       hControl    xfigure object of type listbox or dropdown
%       pos         numeric array that represents the index of the
%                   requested, replaced, or inserted rows
%       setstr      if given will either replace or insert strings,
%                   for single lines can be a char array, otherwise
%                   a cell array of chars
%       insert      if true inserts strings passed by setstr, otherwise
%                   replace strings
%
% Note: (1) if pos is missing, but setstr is given, setstr will be used to
%           either replace the GUI selected subset or, if insert is true,
%           inserted at the end of the control's string list
%       (2) if pos is given but no setstr is passed, mstring returns
%           the subset of strings indexed by pos
%       (3) without any arguments, mstring returns the GUI selected
%           entries, equivalent to mstring(hControl, hControl.Value)
%       (4) if insert is true, pos gives the lines, where setstr entries
%           are inserted BEFORE, and -Inf and Inf marks the beginning and
%           end of the controls string list
%       (5) whenever a string is inserted, the first inserted entry will
%           be selected after the call is complete

% Version:  v1.1
% Build:    16041811
% Date:     Apr-18 2016, 11:52 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2016, Jochen Weber
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
if numel(xo) ~= 1
    error( ...
        'xfigure:InvalidObjectSize', ...
        'mstring only works for singular objects.' ...
    );
end

% try retrieving/setting/inserting strings
try
    mstr = mstring(xo, varargin{:});
catch xfigerror
    rethrow(xfigerror);
end
