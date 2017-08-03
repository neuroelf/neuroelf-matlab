function icbmout = tal2icbm(talin, rounded)
% tal2icbm  - converts coordinates from TAL brain to segmented ICBM space
%
% FORMAT:       mniout = tal2icbm(talin [, rounded])
%
% Input fields:
%
%       talin       N-by-3 or 3-by-N matrix of coordinates
%       rounded     1x1 double, if given, coordinates are rounded to
%                   specified number of digits
%
% Output fields:
%
%       icbmout     is the coordinate matrix with MNI-space points

% Version:  v0.9b
% Build:    11050712
% Date:     Nov-20 2010, 11:09 PM EST
% Editor:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% persistent transformation matrices
persistent t2i_trf;
if isempty(t2i_trf) || ...
   ~isstruct(t2i_trf)
    try
        t2i_trf = load(neuroelf_file('t', 'talairach_seg_inv_sn.mat'));
        if ~isfield(t2i_trf, 'Affine') || ...
           ~isfield(t2i_trf, 'VG') || ...
           ~isstruct(t2i_trf.VG) || ...
            isempty(t2i_trf.VG) || ...
           ~isfield(t2i_trf.VG, 'dim') || ...
           ~isfield(t2i_trf.VG, 'mat')
            error( ...
                'neuroelf:BadFileContent', ...
                'Talairach->ICBM normalization file invalid or missing.' ...
            );
        end
    catch ne_eo;
        t2i_trf = [];
        rethrow(ne_eo);
    end
    t2i_trf.VG = t2i_trf.VG(1);
    t2i_trf.VGmat = inv(t2i_trf.VG.mat);
end

% argument check
if nargin < 1 || ...
    length(size(talin)) > 2 || ...
   ~isa(talin, 'double')
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing input argument mniin.' ...
    );
end

% transpose ?
dimdim = find(size(talin) == 3);
if isempty(dimdim)
    error( ...
        'neuroelf:BadArguments', ...
        'talin argument must be a N-by-3 or 3-by-N matrix' ...
    );
end

% transpose as needed
if dimdim(1) == 1 && ...
    numel(dimdim) == 1
    talin = talin';
end

% use applyspmsnc for transformation
icbmout = applyspmsnc(talin, t2i_trf.Tr, t2i_trf.VG.dim, ...
    t2i_trf.VGmat, t2i_trf.VF(1).mat * t2i_trf.Affine);

% retranspose to match input
if dimdim(1) == 1 && ...
    numel(dimdim) == 1
    icbmout = icbmout';
end

% round output
if nargin > 1 && ...
    isa(rounded, 'double') && ...
   ~isempty(rounded)
    icbmout = (1/10^rounded(1)) * round(icbmout * 10^rounded(1));
end
