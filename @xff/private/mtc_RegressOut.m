function xo = mtc_RegressOut(xo, varargin)
% MTC::RegressOut  - regress out nuisance variance
%
% FORMAT:       [mtc = ] mtc.RegressOut(sdm [, sdm2, ...] [, opts])
%
% Input fields:
%
%       sdm         design or TxR regressors to regress out
%       sdm2, ...   additional nuisance (multiple files)
%       opts        options settings
%        .tfiltfrq  number of temporal filtering frequencies (default: 0)
%        .tfilttyp  filtering type, either of 'DCT', {'Fourier'}
%        .trans     perform either 'psc' or 'z' transform (default: 'none')
%
% Output fields:
%
%       mtc         MTC with variance regressed out (mean intact)
%
% Using: calcbetas, psctrans, tempfilter, ztrans.

% Version:  v1.1
% Build:    16032210
% Date:     Mar-22 2016, 10:54 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'mtc')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
mtcd = double(bc.MTCData);
mtcs = size(mtcd);
nvol = mtcs(1);
if numel(varargin{end}) ~= 1 || ~isstruct(varargin{end})
    opts = struct;
else
    opts = varargin{end};
end
if ~isfield(opts, 'tfiltfrq') || numel(opts.tfiltfrq) ~= 1 || ~isa(opts.tfiltfrq, 'double') || ...
    isinf(opts.tfiltfrq) || opts.tfiltfrq < 0 || opts.tfiltfrq > 12
    opts.tfiltfrq = 0;
end
if ~isfield(opts, 'tfilttyp') || ~ischar(opts.tfilttyp) || ...
   ~any(strcmpi(opts.tfilttyp(:)', {'dct', 'fourier'}))
    opts.tfilttyp = 'fourier';
else
    opts.tfilttyp = lower(opts.tfilttyp(:)');
end
if ~isfield(opts, 'trans') || ~ischar(opts.trans) || isempty(opts.trans) || ...
   ~any(lower(opts.trans(1)) == 'npz')
    opts.trans = 'n';
else
    opts.trans = lower(opts.trans(1));
end

% create design matrix for regression
X = zeros(nvol, 0);

% add regressors
for ac = 1:nargin-1

    % numeric matrix
    if isnumeric(varargin{ac}) && size(varargin{ac}, 1) == nvol

        % add to X
        X = [X, varargin{ac}(:, :)];

    % RTC/SDM
    elseif numel(varargin{ac}) == 1 && xffisobject(varargin{ac}, true, 'sdm')

        % get content
        sdmc = varargin{ac}.C;

        % size check
        if size(sdmc.SDMMatrix, 1) ~= nvol
            error('neuroelf:xff:badArgument', 'SDM size mismatch with MTCData.');
        end

        % add to design
        X = [X, sdmc.SDMMatrix];

    % text/mat file
    elseif ischar(varargin{ac}) && ~isempty(varargin{ac})

        % try to load
        try
            xa = load(varargin{ac}(:)');
        catch xfferror
            neuroelf_lasterr(xfferror);
            try
                xa = xff(varargin{ac}(:)');
                if xffisobject(xa, true)
                    xac = xa.C;
                    delete(xa);
                else
                    xac =[];
                end
                if isstruct(xac) && isfield(xac, 'SDMMatrix')
                    xa = xac.SDMMatrix;
                    if size(xa, 1) ~= nvol
                        error('neuroelf:xff:badArgument', 'SDM size mismatch with MTCData.');
                    end
                else
                    xa = [];
                end
            catch xfferror
                neuroelf_lasterr(xfferror);
                xa = [];
            end
        end
        if isnumeric(xa) && size(xa, 1) == nvol
            X = [X, xa(:, :)];
        elseif isstruct(xa) && numel(fieldnames(xa)) == 1
            xf = fieldnames(xa);
            xa = xa.(xf{1});
            if isnumeric(xa) && size(xa, 1) == nvol
                X = [X, xa(:, :)];
            end
        end
    end
end

% filter content
if opts.tfiltfrq > 0

    % prepare tempfilter options
    topts = opts;
    topts.spat = false;
    topts.tdim = 1;
    topts.temp = true;
    topts.tempdt = false;
    if opts.tfilttyp(1) == 'd'
        if opts.tfiltfrq > 0
            topts.tempdct = ceil(0.001 * bc.TR * nvol / opts.tfiltfrq);
        else
            topts.tempdct = Inf;
        end
        topts.tempsc = 0;
    else
        topts.tempdct = Inf;
        topts.tempsc = opts.tfiltfrq;
    end

    % temp filter data of first object
    [null, Xf] = ne_methods.tempfilter(zeros(nvol, 1), topts);

    % and temp filter regressors
    X = [X, Xf];
end

% cleanup
X(:, any(isnan(X) | isinf(X))) = [];
X(:, sum(abs(diff(X))) < sqrt(eps)) = [];

% constant
X(:, end + 1) = 1;

% calcbetas
b = ne_methods.calcbetas(X, mtcd);

% set constant to 0 (unless z-trans)
if opts.trans ~= 'z'
    b(:, end) = 0;
end

% regress out
mtcd = mtcd - X * b';

% transformation
if opts.trans == 'p'
    mtcd = ne_methods.psctrans(mtcd);
elseif opts.trans == 'z'
    mtcd = ne_methods.ztrans(mtcd);
end

% remove infs/nans
mtcd(repmat(any(isinf(mtcd) | isnan(mtcd)), nvol, 1)) = 0;

% set content
bc.MTCData = single(mtcd);
xo.C = bc;
