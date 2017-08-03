function xo = vmr_Update(xo, F, S, V)
% VMR::Update  - called after subsasgn for VMRs
%
% Using: makelabel.

% Version:  v1.1
% Build:    16021413
% Date:     Feb-14 2016, 1:04 PM EST
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
makelabel = ne_methods.makelabel;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmr') || ~ischar(F) || isempty(F)
    error('neuroelf:xff:BadArgument', 'Invalid call to %s.', mfilename);
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
F = makelabel(F(:)');

% F valid
if ~isfield(bc, F)
    error('neuroelf:xff:invalidProperty', 'Cannot find property ''%s'' for type VMR.', F);
end

% what field has changed
switch (lower(F))

    % reframing
    case {'dimx', 'dimy', 'dimz'}

        % which target type
        if bc.VMR8bit
            is8bit = true;
        else
            is8bit = false;
        end

        % get destination and current data size
        dsize = [bc.DimX, bc.DimY, bc.DimZ];
        csize = size(bc.VMRData);

        % changes ?
        if any(dsize ~= csize)

            % dim1
            if dsize(1) < csize(1)
                srange1 = round(0.5 + (csize(1) - dsize(1)) / 2);
                srange1(2) = srange1 + dsize(1) - 1;
                trange1 = [1, dsize(1)];
            elseif dsize(1) > csize(1)
                srange1 = [1, csize(1)];
                trange1 = round(0.5 + (dsize(1) - csize(1)) / 2);
                trange1(2) = trange1 + csize(1) - 1;
            else
                srange1 = [1, csize(1)];
                trange1 = srange1;
            end

            % dim2
            if dsize(2) < csize(2)
                srange2 = round(0.5 + (csize(2) - dsize(2)) / 2);
                srange2(2) = srange2 + dsize(2) - 1;
                trange2 = [1, dsize(2)];
            elseif dsize(2) > csize(2)
                srange2 = [1,csize(2)];
                trange2 = round(0.5 + (dsize(2) - csize(2)) / 2);
                trange2(2) = trange2 + csize(2) - 1;
            else
                srange2 = [1, csize(2)];
                trange2 = srange2;
            end

            % dim3
            if dsize(3) < csize(3)
                srange3 = round(0.5 + (csize(3) - dsize(3)) / 2);
                srange3(2) = srange3 + dsize(3) - 1;
                trange3 = [1, dsize(3)];
            elseif dsize(3) > csize(3)
                srange3 = [1,csize(3)];
                trange3 = round(0.5 + (dsize(3) - csize(3)) / 2);
                trange3(2) = trange3 + csize(3) - 1;
            else
                srange3 = [1, csize(3)];
                trange3 = srange3;
            end

            % get source data
            VMRsrc = bc.VMRData( ...
                srange1(1):srange1(2), ...
                srange2(1):srange2(2), ...
                srange3(1):srange3(2));

            % create target data
            if is8bit
                bc.VMRData = uint8(0);
            else
                bc.VMRData = uint16(0);
            end
            bc.VMRData(1:dsize(1), 1:dsize(2), 1:dsize(3)) = 0;

            % put source in target
            bc.VMRData( ...
                trange1(1):trange1(2), ...
                trange2(1):trange2(2), ...
                trange3(1):trange3(2)) = VMRsrc;
        end


    % changing file version
    case 'fileversion'

        % get wanted version
        reqv = bc.FileVersion;

        % check version
        if isempty(reqv) || ...
           ~isa(reqv, 'double') || ...
            isnan(reqv(1)) || ...
            isinf(reqv(1)) || ...
            fix(reqv(1)) ~= reqv(1) || ...
            reqv(1) < 0 || ...
            reqv(1) > 3

            % give warning if required
            warning('neuroelf:xff:invalidPropertyValue', 'Invalid FileVersion value given.');

            % set better value and go on
            reqv = 3;
            bc.FileVersion = reqv;
        end

    % set to 8 bit (max intensity 223)
    case 'vmr8bit'

        % get target state
        tt = bc.VMR8bit;
        if isempty(tt)
            tt = false;
        elseif ~isnumeric(tt)
            tt = true;
        elseif tt(1)
            tt = true;
        else
            tt = false;
        end

        % only check type for 8bit (other is checked in VMRData handler)
        if tt && ~isa(bc.VMRData, 'uint8')

            % set and find maxval
            maxval = 223;
            curmax = double(max(bc.VMRData(:)));

            % only reset if greater than possible values!
            if curmax > 255

                % recalculate intensities
                bc.VMRData = uint8(round(double(bc.VMRData) * (maxval/curmax)));

            % else just adapt datatype
            elseif ~isa(bc.VMRData, 'uint8')

                % set datatype
                bc.VMRData = uint8(bc.VMRData);
            end
        end

        % if ~tt and ~uint16
        if ~tt && ~isa(bc.VMRData, 'uint16')

            % set datatype
            bc.VMRData = uint16(bc.VMRData);
        end


    % new data set (either partial or complete!)
    case 'vmrdata'

        % get current datatype
        if bc.VMR8bit
            is8bit = true;
            maxval = uint8(255);
        else
            is8bit = false;
            maxval = uint16(4095);
        end

        % check type
        checkcont = false;
        if is8bit
            if ~isa(bc.VMRData, 'uint8')
                checkcont = true;
            end
        else
            if ~isa(bc.VMRData, 'uint16')
                checkcont = true;
            end
        end

        % if type must be converted
        if checkcont

            % check numeric content
            if ~isnumeric(bc.VMRData)

                % give warning
                warning('neuroelf:xff:invalidDataType','Invalid datatype for VMRData property.');

                % and fill with type according zeros
                csize = [bc.DimX, bc.DimY, bc.DimZ];
                if is8bit
                    bc.VMRData = uint8(0);
                else
                    bc.VMRData = uint16(0);
                end
                bc.VMRData(1:csize(1), 1:csize(2), 1:csize(3)) = 0;
            else

                % replace NaN's and Inf's
                bc.VMRData(isnan(bc.VMRData)) = 0;
                bc.VMRData(isinf(bc.VMRData)) = maxval;

                % get min and max val
                curmin = double(min(bc.VMRData(:)));
                if curmin > 0
                    curmin = 0;
                end
                curmax = double(max(bc.VMRData(:)));

                % reframe intensities
                if is8bit
                    if curmax > 255
                        maxval = uint8(223);
                    end
                    bc.VMRData = uint8(round( ...
                        (double(bc.VMRData) - curmin) / ...
                        (curmax + eps - curmin) * double(maxval)));
                else
                    bc.VMRData = uint16(round( ...
                        (double(bc.VMRData) - curmin) / ...
                        (curmax + eps - curmin) * double(maxval)));
                end
            end
        end

        % get source and current data size
        ssize = [bc.DimX, bc.DimY, bc.DimZ];
        csize = size(bc.VMRData);

        % changes ?
        if any(ssize ~= csize)
            bc.DimX = csize(1);
            bc.DimY = csize(2);
            bc.DimZ = csize(3);
        end
end

% update contents
xo.C = bc;
