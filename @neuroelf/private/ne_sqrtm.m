function w = ne_sqrtm(v)
% ne_sqrtm  - matrix sqrt for symmetric, positive semi-definite matrices
%
% FORMAT:       w = ne_sqrtm(v)
%
% Input fields:
%
%       v           symmetric, positive semi-definite matrix
%
% Output fields:
%
%       w           (weighting) matrix fullfilling: w' * w = inv(v)
%
% This routine covers and extends sqrtm functionality by using a
% computationally expedient approximation that can handle sparse
% symmetric positive semi-definite matrices.

% Version:  v0.0
% Build:    13020214
% Date:     Feb-02 2013, 2:49 PM EST
% Author:   Karl Friston
% Editor:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_sqrtm.m 4124 2010-11-18 16:56:53Z karl $

%--------------------------------------------------------------------------
if issparse(v)
    v = full(v);
end
[u, s] = svd(inv(v), 0);
w = u * diag(sqrt(diag(s))) * u';
