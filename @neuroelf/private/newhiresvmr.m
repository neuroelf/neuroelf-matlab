function hvmr = newhiresvmr(d)
% newhiresvmr  - create new hires (0.5mm) VMR object
%
% FORMAT:       hvmr = newhiresvmr([d])
%
% Input fields:
%
%       d           if given and 512x512x512 uint8, put into data
%
% Output fields:
%
%       hvmr        0.5mm resolution VMR (only 8-bit VMRData)

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:29 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2016, Jochen Weber
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

% create VMR
hvmr = xff('new:vmr');

% settings
xffroot = xff();
updstate = xffroot.UpdateState('vmr', false);

% fill VMR
hvmr.VMRData = uint8(0);
hvmr.FramingCube = 512;
hvmr.PosInfoVerified = 1;
hvmr.Slice1CenterX = -127.75;
hvmr.SliceNCenterX = 127.75;
hvmr.NRows = 512;
hvmr.NCols = 512;
hvmr.SliceThickness = 0.5;
hvmr.VoxResX = 0.5;
hvmr.VoxResY = 0.5;
hvmr.VoxResZ = 0.5;
hvmr.VoxResInTalairach = 1;
hvmr.VoxResVerified = 1;
hvmr.VMRData16 = uint16([]);
if nargin > 0 && ...
    ndims(d) == 3 && ...
    isa(d, 'uint8') && ...
   ~any(size(d) > 512)
    hvmr.VMRData = d;
else
    hvmr.VMRData(512, 512, 512) = uint8(0);
end

% reset updatestate flag
xffroot.UpdateState('vmr', updstate);
