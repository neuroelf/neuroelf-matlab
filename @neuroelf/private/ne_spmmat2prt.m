% FUNCTION ne_spmmat2prt: perform SPM.mat -> PRT conversion
function varargout = ne_spmmat2prt(varargin)

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-10 2011, 4:52 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% request SPM.mat filename
[spmmatfile, spmmatpath] = uigetfile({'SPM*.mat', 'SPM design matrix files'}, ...
    'Please select SPM.mat file to convert to PRT...', 'SPM*.mat', 'MultiSelect', 'off');
if isequal(spmmatfile, 0) || ...
    isequal(spmmatpath, 0) || ...
    exist([spmmatpath '/' spmmatfile], 'file') ~= 2
    return;
end

% request protocol target file(s)
[prtfile, prtpath] = uiputfile({'*.prt', 'BrainVoyager Protocol files'}, ...
    'Please specify the first target filename...', 'SPM*.prt');
if isequal(prtfile, 0) || ...
    isequal(prtpath, 0)
    return;
end

% conversion
try
    spmmatfile = strrep(strrep([spmmatpath '/' spmmatfile], '\', '/'), '//', '/');
    prtfile = strrep(strrep([prtpath '/' prtfile], '\', '/'), '//', '/');
    [p, ps] = spmmat2prt(spmmatfile, prtfile);
catch ne_eo;
    errordlg(['Error converting SPM.mat -> PRT: ' ne_eo.message], 'NeuroElf GUI - error', 'modal');
    return;
end

% clear protocol objects
clearxffobjects(ps);
