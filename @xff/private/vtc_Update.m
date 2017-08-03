function xo = vtc_Update(xo, F, S, V)
% VTC::Update  - called after subsasgn for VTCs
%
% Using: makelabel, minmaxmean.

% Version:  v1.1
% Build:    16021320
% Date:     Feb-13 2016, 8:35 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2016, Jochen Weber
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
minmaxmean = ne_methods.minmaxmean;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vtc') || ~ischar(F) || isempty(F)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 3 || ~isstruct(S)
    S = struct;
    S.type = '.';
    S.subs = F;
end
if nargin < 4
    V = [];
end

% get content
bc = xo.C;

% linearize
F = ne_methods.makelabel(F(:)');

% F valid
if ~isfield(bc, F)
    error('neuroelf:xff:invalidProperty', 'Cannot find property ''%s'' for type VTC.', F);
end

% what field has changed
switch (lower(F))

    % DataType
    case 'datatype'

        % get requested DataType
        dt = bc.DataType;
        switch (dt)

            % the uint16 DataType
            case 1

                % check for transio
                if istransio(bc.VTCData)

                    % check class of transio
                    dv = bc.VTCData(1);

                    % if not uint16, transform into data matrix in mem
                    if ~isa(dv, 'uint16')

                        % get data as single
                        dv = single(resolve(bc.VTCData));

                        % get min max mean
                        mmm = minmaxmean(dv(:));

                        % remove min
                        dv = dv - mmm(1);
                        mmm = mmm - mmm(1);

                        % check if up-scaling is a good idea
                        if mmm(2) <= 16383.5 && mmm(2) > 0

                            % find best integer factor
                            sf = floor(32767 / mmm(2));

                            % up-scale data
                            dv = sf .* dv;

                        % down-scaling instead ?
                        elseif mmm(2) > 32767

                            % just do it
                            dv = (32767 / mmm(2)) .* dv;
                        end

                        % round and uint16 data in mem
                        bc.VTCData = uint16(round(dv));
                    end

                % if the data is not transio and not uint 16
                elseif ~isa(bc.VTCData, 'uint16')

                    % for int8
                    switch (class(bc.VTCData))
                        case {'int16'}
                            mmm = minmaxmean(bc.VTCData(:));
                            bc.VTCData = uint16(single(bc.VTCData) - mmm(1));
                        case {'int8'}
                            bc.VTCData = uint16(128 + single(bc.VTCData));
                        case {'logical', 'uint8'}
                            bc.VTCData = uint16(bc.VTCData);
                        otherwise

                            % follow transio logic, get data as single
                            dv = double(bc.VTCData);

                            % get min max mean
                            mmm = minmaxmean(dv(:));

                            % remove min
                            dv = dv - mmm(1);
                            mmm = mmm - mmm(1);

                            % check if up-scaling is a good idea
                            if mmm(2) <= 16383.5 && ...
                                mmm(2) > 0

                                % find best integer factor
                                sf = floor(32767 / mmm(2));

                                % up-scale data
                                dv = sf .* dv;

                            % down-scaling instead ?
                            elseif mmm(2) > 32767

                                % just do it
                                dv = (32767 / mmm(2)) .* dv;
                            end

                            % round and uint16 data in mem
                            bc.VTCData = uint16(round(dv));
                    end
                end

            case 2

                % also set FileVersion if < 3
                if bc.FileVersion < 3
                    bc.FileVersion = 3;
                end

                % and make sure data is single
                if istransio(bc.VTCData)
                    dv = bc.VTCData(1);
                    if ~isa(dv, 'single')
                        bc.VTCData = single(resolve(bc.VTCData));
                    end
                elseif ~isa(bc.VTCData, 'single')
                    bc.VTCData = single(bc.VTCData);
                end

            otherwise
                warning('neuroelf:xff:invalidSetting', ...
                    'Invalid DataType setting made. Reverting...');
        end
end

% update contents
xo.C = bc;
