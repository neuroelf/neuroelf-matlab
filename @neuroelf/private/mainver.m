function mv = mainver
% mainver  - get matlab main version
%
% FORMAT:       mv = mainver
%
% Output Fields:
%
%       mv          numeric main version number
%
% See also ver, version

% Version:  v0.9a
% Build:    11050711
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

% get version info?
persistent mainver_ver;
if isempty(mainver_ver)
    mv_temp = ver('matlab');
    if isempty(mv_temp) || ...
       ~isstruct(mv_temp)
        mv_temp = ver('octave');
        if isstruct(mv_temp) && ...
            numel(mv_temp) == 1
            mainver_ver = str2double(mv_temp(1).Version(1));
        else
            warning( ...
                'neuroelf:InternalError', ...
                'Error getting MATLAB version number.' ...
            );
            mainver_ver = 6;
        end
        mv = mainver_ver;
        return;
    end
    mv_temp = mv_temp(1).Version;
    mv_ivchar = find(mv_temp < 48 | mv_temp > 57);
    if ~isempty(mv_ivchar)
        mv_temp(mv_ivchar(1):end) = [];
    end
    if ~isempty(mv_temp)
        mainver_ver = str2double(mv_temp);
    else
        mainver_ver = 6;
        warning( ...
            'neuroelf:InternalError', ...
            'Error getting MATLAB version number.' ...
        );
    end
end

% return version info
mv = mainver_ver;
