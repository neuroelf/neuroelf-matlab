function [varargout] = ne_fastica(mixedsig, opts)
% ne_fastica  - perform ICA with FastICA algorithm
%
% FORMAT:       icasig = ne_fastica(mixedsig [, ...])
%    -or-       [icasig, A, W] = ne_fastica(mixedsig [, ...])
%    -or-       [A, W] = ne_fastica(mixedsig [, ...])
%
% Input fields:
%
%       mixedsig    GxS double G signals by S samples mixed signals
%                   (each row is one mixed signal)
%       opts        struct with additional settings
%        .a1        parameter (number), used when g := 'tanh', default: 1
%        .a2        parameter (number), used when g := 'gaus', default: 1
%        .approach  decorrelation approach used, either 'symm' or {'defl'}
%                   whereas symmetric (symm) estimates all of the
%                   independent components in parallel and deflation (defl)
%                   estimates the independent components one-by-one, like
%                   in projection pursuit; default is defl
%        .eig1      first (largest) eigenvalue to be retained (default: 1)
%        .eign      last (smallest) eigenvalue to be retained (default: all)
%        .epsilon   stopping criterion (number, default: 0.0001)
%        .g         nonlinearity g in the fixed-point algorithm, one of
%                   {'pow3'}    g(u) = u .^ 3
%                    'tanh'     g(u) = tanh(a1 * u)
%                    'gauss'    g(u) = u * exp(-a2 * u .^ 2 ./ 2)
%                    'skew'     g(u) = u .^ 2
%        .gfinetune g when fine-tuning, same selection as g, plus {'off'}
%        .iguess    initial guess for A (default: random)
%                   this allows for a "one-more" approach:
%                   [ica, A, W] = ne_fastica(mix, struct('numics',3));
%                   [ica2, A2, W2] = ne_fastica(mix, struct('iguess', A, 'numics', 4));
%        .maxftune  max number of iterations in fine-tuning (default: 100)
%        .maxiter   maximum number of iterations (number, default: 1000)
%        .mu        step size (number), default: 1
%                   if the value of mu is not 1 the program will use the
%                   stabilized version of the algorithm (see also .stab)
%        .numics    number of components to estimate (default: data dim)
%        .smpsize   percentage of samples used ([0 .. 1], default: 1)
%        .stab      boolean flag, use stabilized algorithm (default: false)
%                   if stabilization is enabled the value of mu can be
%                   temporarily be halved in case the program detects that
%                   the algorithm is stuck between two points (so-called
%                   stroke); also if there is no convergence before half
%                   of the maximum number of iterations has been reached
%                   then mu will be halved for the rest of the rounds
%        .step      if set, either of 'pca' or 'white', compute only
%                   PCA part (possibly with whitening of data)
%                   examples are:
%                   [e, d] = ne_fastica(mixedsig, struct('step', 'pca'));
%                   [ws, WM, DWM] = ne_fastica(mixedsig, struct('step', white'));
%
% Output fields:
%
%       icasig      rows contain the estimated independent components
%       A           mixing matrix A
%       W           separating matrix W
%       e, d        eigenvalues and eigenvectors (for step := 'pca')
%
% Note: this algorithm (code) was re-used by permission from the FastICA
%       Matlab package version 2.5 to be found at
%       http://research.ics.tkk.fi/ica/fastica/
%       and is Copyright (c) by
%       Hugo Gävert, Jarmo Hurri, Jaakko Särelä, and Aapo Hyvärinen.
%
% Note: the function name was deliberately altered so as not to clash!

% Examples:
%
%       icasig = ne_fastica(mixedsig, struct('approach', 'symm', 'g', 'tanh'));
%           do ICA with tanh nonlinearity and in parallel (like
%           maximum likelihood estimation for supergaussian data).
%
%       icasig = ne_fastica(mixedsig, struct('eign', 10, 'numics', 3));
%           first, reduce dimension to 10, and then estimate only 3
%           independent components.
%
% @(#)$Id: fastica.m,v 1.14 2005/10/19 13:05:34 jarmo Exp $ (derivate)

% Version:  v1.0
% Build:    15050511
% Date:     May-05 2015, 11:10 AM EST
% Author:   Hugo Gävert, Jarmo Hurri, Jaakko Särelä, and Aapo Hyvärinen
% Editor:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2015, Jochen Weber

% argument check
if nargin < 1 || ...
   ~isnumeric(mixedsig) || ...
    ndims(mixedsig) > 2
    error( ...
        'neuroelf:FastICAError', ...
        'You must supply the mixed data as input argument.' ...
    );
end
mixedsig = double(mixedsig);
if any(isinf(mixedsig(:)) | isnan(mixedsig(:)))
    error( ...
        'neuroelf:FastICAError', ...
        'The mixed data must not contain Inf/NaN values.' ...
    );
end
[nsig, nsmp] = size(mixedsig);
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'a1') || ...
   ~isa(opts.a1, 'double') || ...
    numel(opts.a1) ~= 1 || ...
    isinf(opts.a1) || ...
    isnan(opts.a1) || ...
    opts.a1 <= 0
    opts.a1 = 1;
end
if ~isfield(opts, 'a2') || ...
   ~isa(opts.a2, 'double') || ...
    numel(opts.a2) ~= 1 || ...
    isinf(opts.a2) || ...
    isnan(opts.a2) || ...
    opts.a2 <= 0
    opts.a2 = 1;
end
if ~isfield(opts, 'approach') || ...
   ~ischar(opts.approach) || ...
    isempty(opts.approach) || ...
   ~any(strcmpi(opts.approach(:)', ...
        {'d', 'defl', 'deflate', 's', 'symm', 'symmetric'}))
    opts.approach = 2;
else
    switch (lower(opts.approach(1)))
        case {'d'}
            opts.approach = 2;
        otherwise
            opts.approach = 1;
    end
end
if ~isfield(opts, 'eig1') || ...
   ~isa(opts.eig1, 'double') || ...
    numel(opts.eig1) ~= 1 || ...
    isinf(opts.eig1) || ...
    isnan(opts.eig1) || ...
    opts.eig1 < 1 || ...
    opts.eig1 >= nsig || ...
    opts.eig1 ~= fix(opts.eig1)
    opts.eig1 = 1;
end
if ~isfield(opts, 'eign') || ...
   ~isa(opts.eign, 'double') || ...
    numel(opts.eign) ~= 1 || ...
    isinf(opts.eign) || ...
    isnan(opts.eign) || ...
    opts.eign < opts.eig1 || ...
    opts.eign >= nsig || ...
    opts.eign ~= fix(opts.eign)
    opts.eign = nsig;
end
if ~isfield(opts, 'epsilon') || ...
   ~isa(opts.epsilon, 'double') || ...
    numel(opts.epsilon) ~= 1 || ...
    isinf(opts.epsilon) || ...
    isnan(opts.epsilon) || ...
    opts.epsilon <= 0
    opts.epsilon = 0.0001;
end
if ~isfield(opts, 'g') || ...
   ~ischar(opts.g) || ...
    isempty(opts.g) || ...
   ~any(strcmpi(opts.g(:)', {'gauss', 'pow3', 'skew', 'tanh'}))
    opts.g = 'pow3';
else
    opts.g = lower(opts.g(:)');
end
if ~isfield(opts, 'gfinetune') || ...
   ~ischar(opts.gfinetune) || ...
    isempty(opts.gfinetune) || ...
   ~any(strcmpi(opts.gfinetune(:)', {'gauss', 'off', 'pow3', 'skew', 'tanh'}))
    opts.gfinetune = 'off';
else
    opts.gfinetune = lower(opts.gfinetune(:)');
end
if ~isfield(opts, 'iguess') || ...
   ~isa(opts.iguess, 'double') || ...
    ndims(opts.iguess) > 2 || ...
   ~all(size(opts.iguess) == nsig) || ...
    any(isinf(opts.iguess(:)) | isnan(opts.iguess(:)))
    opts.iguess = 0;
    opts.istate = false;
else
    opts.istate = true;
end
if ~isfield(opts, 'maxftune') || ...
   ~isa(opts.maxftune, 'double') || ...
    numel(opts.maxftune) ~= 1 || ...
    isinf(opts.maxftune) || ...
    isnan(opts.maxftune) || ...
    opts.maxftune < 0
    opts.maxftune = 5;
end
if ~isfield(opts, 'maxiter') || ...
   ~isa(opts.maxiter, 'double') || ...
    numel(opts.maxiter) ~= 1 || ...
    isinf(opts.maxiter) || ...
    isnan(opts.maxiter) || ...
    opts.maxiter <= 0
    opts.maxiter = 1000;
else
    opts.maxiter = fix(opts.maxiter);
end
if ~isfield(opts, 'mu') || ...
   ~isa(opts.mu, 'double') || ...
    numel(opts.mu) ~= 1 || ...
    isinf(opts.mu) || ...
    isnan(opts.mu) || ...
    opts.mu <= 0
    opts.mu = 1;
end
if ~isfield(opts, 'numics') || ...
   ~isa(opts.numics, 'double') || ...
    numel(opts.numics) ~= 1 || ...
    isinf(opts.numics) || ...
    isnan(opts.numics) || ...
    opts.numics <= 0 || ...
    opts.numics > nsig
    opts.numics = nsig;
else
    opts.numics = fix(opts.numics);
end
if ~isfield(opts, 'pbar') || ...
    numel(opts.pbar) ~=1 || ...
   ~any(strcmpi(class(opts.pbar), {'xfigure', 'xprogress'}))
    opts.pbar = [];
end
if ~isfield(opts, 'smpsize') || ...
   ~isa(opts.smpsize, 'double') || ...
    numel(opts.smpsize) ~= 1 || ...
    isinf(opts.smpsize) || ...
    isnan(opts.smpsize)
    opts.smpsize = 1;
else
    opts.smpsize = min(1, max(2 / nsig, opts.smpsize));
end
if ~isfield(opts, 'stab') || ...
   ~islogical(opts.stab) || ...
    numel(opts.stab) ~= 1
    opts.stab = false;
end
if ~isfield(opts, 'step') || ...
   ~ischar(opts.step) || ...
    isempty(opts.step) || ...
   ~any(strcmpi(opts.step(:)', ...
        {'a', 'all', 'f', 'full', 'p', 'pca', 'w', 'white'}))
    opts.step = 3;
else
    switch (lower(opts.step(1)))
        case {'p'}
            opts.step = 1;
        case {'w'}
            opts.step = 2;
        otherwise
            opts.step = 3;
    end
end

% progress par?
pbar = [];
if isempty(opts.pbar) && ...
    opts.step > 2
    try
        pbar = xprogress;
        xprogress(pbar, 'setposition', [80, 200, 640, 36]);
        xprogress(pbar, 'settitle', ...
            sprintf('Performing ICA on %dx%d matrix...', size(mixedsig)));
        xprogress(pbar, 0, 'Removing mean...', 'visible', 0, 1);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        pbar = [];
    end
elseif opts.step > 2
    pbar = opts.pbar;
    pbar.Progress(pbmin, 'Removing mean...');
end

% remove mean from signal
[mixedsig, mixedmean] = fi_remmean(mixedsig);

% calculate PCA
if ~isempty(pbar)
    pbar.Progress(0.1 / min(size(mixedsig)), 'Performing PCA...');
end
try
    [E, D] = fi_pcamat(mixedsig, opts.eig1, opts.eign);
catch ne_eo;
    if ~isempty(pbar)
        closebar(pbar)
    end
    rethrow(ne_eo);
end

% only PCA
if opts.step < 2

    % return PCA results
    varargout{1} = E;
    varargout{2} = D;
    return;
end

% whitening the data
if ~isempty(pbar)
    pbar.Progress(0.25 / min(size(mixedsig)), 'Whitening the data...');
end
[X, whiteningMatrix, dewhiteningMatrix] = fi_whitenv(mixedsig, E, D);
if ~isreal(X)
    error( ...
        'neuroelf:FastICAError', ...
        'Input has an imaginary part.' ...
    );
end

% no further processing
if opts.step < 3

    % what output
    if nargout == 2
        varargout{1} = whiteningMatrix;
        varargout{2} = dewhiteningMatrix;
    else
        varargout{1} = X;
        if nargout > 1
            varargout{2} = whiteningMatrix;
            varargout{3} = dewhiteningMatrix;
        end
    end
    return;
end

% the number of IC's must be less or equal to the dimension of data
nsig = size(X, 1);
if opts.numics > nsig
    opts.numics = nsig;
end
pstep = 1 / (opts.numics + 1);

% now comes the main portion of the code...

% checking the value for numics
[numvecs, numwsmp] = size(X);
if numvecs < opts.numics
    error( ...
        'neuroelf:FastICAError', ...
        'Whitened signal must have dim >= numics.' ...
    );
end

% better range for sample size
if opts.smpsize < 1 && ...
   (opts.smpsize * numwsmp) < 1000
    opts.smpsize = min(1000 / numwsmp, 1);
end

% Checking the value for nonlinearity.
switch (opts.g)
    case {'pow3'}
        gOrig = 10;
    case {'tanh'}
        gOrig = 20;
    case {'gauss'}
        gOrig = 30;
    case {'skew'}
        gOrig = 40;
end
if opts.smpsize ~= 1
    gOrig = gOrig + 2;
end
if opts.mu ~= 1
    gOrig = gOrig + 1;
end
finetune = true;
switch (opts.gfinetune)
    case 'pow3'
        gFine = 10 + 1;
    case 'tanh'
        gFine = 20 + 1;
    case {'gauss'}
        gFine = 30 + 1;
    case 'skew'
        gFine = 40 + 1;
    case 'off'
        if opts.mu ~= 1
            gFine = gOrig;
        else
            gFine = gOrig + 1;
        end
        finetune = false;
end
if ~opts.stab && ...
    opts.mu ~= 1
    opts.stab = true;
end

% some other parameters
myy = opts.mu;
myyOrig = myy;
myyK = 0.01;
failureLimit = 5;

% checking the value for initial state.
if opts.istate
    if size(opts.iguess, 1) ~= size(whiteningMatrix, 2)
        opts.istate = false;
    else
        if size(opts.iguess,2) < opts.numics
            opts.iguess(:, size(opts.iguess, 2) + 1:opts.numics) = ...
                rand(numvecs, opts.numics - size(opts.iguess, 2)) - 0.5;
        elseif size(opts.iguess, 2) > opts.numics
            opts.iguess = opts.iguess(:, 1:opts.numics);
        end
    end
end

% symmetric approach
if opts.approach == 1

    % set some parameters more...
    usedNlinearity = gOrig;
    stroke = 0;
    notFine = true;
    long = false;

    % pre-set A (dewhitened basis vectors)
    A = zeros(numvecs, opts.numics);

    % no initial guess
    if ~opts.istate

        % take random orthonormal initial vectors
        B = orth(randn(numvecs, opts.numics));

    % initial guess given
    else

        % use the given initial vector as the initial state
        B = whiteningMatrix * opts.iguess;
    end

    % some more arrays ...
    BOld = zeros(size(B));
    BOld2 = zeros(size(B));

    % this is the actual fixed-point iteration loop
    solutionfound = false;
    for round = 1:maxiter

        % Symmetric orthogonalization.
        B = B * real(inv(B' * B) ^ 0.5);

        % test for termination condition; note that we consider opposite
        % directions here as well
        minAbsCos = min(abs(diag(B' * BOld)));
        minAbsCos2 = min(abs(diag(B' * BOld2)));

        if (1 - minAbsCos < opts.epsilon)
            if finetune && ...
                notFine
                notFine = 0;
                usedNlinearity = gFine;
                myy = myyK * myyOrig;
                BOld = zeros(size(B));

            else

                % calculate the de-whitened vectors
                A = dewhiteningMatrix * B;
                solutionfound = true;
                break;
            end

        elseif opts.stab
            if stroke == 0 && ...
               (1 - minAbsCos2 < opts.epsilon)
                stroke = myy;
                myy = 0.5 * myy;
                if mod(usedNlinearity, 2) == 0
                    usedNlinearity = usedNlinearity + 1;
                end
            elseif stroke > 0
                myy = stroke;
                stroke = 0;
                if myy == 1 && ...
                    mod(usedNlinearity, 2) ~= 0
                    usedNlinearity = usedNlinearity - 1;
                end
            elseif ~long && ...
                round > (0.5 * maxiter)
                long = true;
                myy = 0.5 * myy;
                if mod(usedNlinearity, 2) == 0
                    usedNlinearity = usedNlinearity + 1;
                end
            end
        end

        % keep history
        BOld2 = BOld;
        BOld = B;

        % sub-approach settings
        switch usedNlinearity

            % pow3
            case {10}
                B = (X * ((X' * B) .^ 3)) / numwsmp - 3 * B;
            case {11}
                % optimoitu - epsilonin kokoisia eroja
                % tämä on optimoitu koodi, katso vanha koodi esim.
                % aikaisemmista versioista kuten 2.0 beta3
                Y = X' * B;
                Gpow3 = Y .^ 3;
                Beta = sum(Y .* Gpow3);
                D = diag(1 ./ (Beta - 3 * numwsmp));
                B = B + myy * B * (Y' * Gpow3 - diag(Beta)) * D;
            case {12}
                Xsub = X(:, fi_getsmp(numwsmp, opts.smpsize));
                B = (Xsub * ((Xsub' * B) .^ 3)) / size(Xsub, 2) - 3 * B;
            case {13}
                % optimoitu
                Ysub = X(:, fi_getsmp(numwsmp, opts.smpsize))' * B;
                Gpow3 = Ysub .^ 3;
                Beta = sum(Ysub .* Gpow3);
                D = diag(1 ./ (Beta - 3 * size(Ysub', 2)));
                B = B + myy * B * (Ysub' * Gpow3 - diag(Beta)) * D;

            % tanh
            case {20}
                hypTan = fi_tanh(a1 * X' * B);
                B = X * hypTan / numwsmp - ones(size(B, 1), 1) * ...
                    sum(1 - hypTan .^ 2) .* B / numwsmp * opts.a1;
            case {21}
                % optimoitu - epsilonin kokoisia
                Y = X' * B;
                hypTan = fi_tanh(opts.a1 * Y);
                Beta = sum(Y .* hypTan);
                D = diag(1 ./ (Beta - opts.a1 * sum(1 - hypTan .^ 2)));
                B = B + myy * B * (Y' * hypTan - diag(Beta)) * D;
            case {22}
                Xsub = X(:, fi_getsmp(numwsmp, opts.smpsize));
                hypTan = fi_tanh(opts.a1 * Xsub' * B);
                B = Xsub * hypTan / size(Xsub, 2) - ones(size(B, 1), 1) * ...
                    sum(1 - hypTan .^ 2) .* B / size(Xsub, 2) * opts.a1;
            case {23}
                % optimoitu
                Y = X(:, fi_getsmp(numwsmp, opts.smpsize))' * B;
                hypTan = fi_tanh(opts.a1 * Y);
                Beta = sum(Y .* hypTan);
                D = diag(1 ./ (Beta - opts.a1 * sum(1 - hypTan .^ 2)));
                B = B + myy * B * (Y' * hypTan - diag(Beta)) * D;

            % gauss
            case {30}
                U = X' * B;
                Usquared = U .^ 2;
                ex = exp(-opts.a2 * Usquared / 2);
                gauss =  U .* ex;
                dGauss = (1 - opts.a2 * Usquared) .* ex;
                B = X * gauss / numwsmp - ones(size(B, 1), 1) * ...
                    sum(dGauss) .* B / numwsmp ;
            case {31}
                % optimoitu
                Y = X' * B;
                ex = exp(-opts.a2 * (Y .^ 2) / 2);
                gauss = Y .* ex;
                Beta = sum(Y .* gauss);
                D = diag(1 ./ (Beta - sum((1 - opts.a2 * (Y .^ 2)) .* ex)));
                B = B + myy * B * (Y' * gauss - diag(Beta)) * D;
            case {32}
                Xsub = X(:, fi_getsmp(numwsmp, opts.smpsize));
                U = Xsub' * B;
                Usquared = U .^ 2;
                ex = exp(-opts.a2 * Usquared / 2);
                gauss =  U .* ex;
                dGauss = (1 - opts.a2 * Usquared) .* ex;
                B = Xsub * gauss / size(Xsub, 2) - ones(size(B, 1), 1) * ...
                    sum(dGauss) .* B / size(Xsub,2);
            case {33}
                % optimoitu
                Y = X(:, fi_getsmp(numwsmp, opts.smpsize))' * B;
                ex = exp(-opts.a2 * (Y .^ 2) / 2);
                gauss = Y .* ex;
                Beta = sum(Y .* gauss);
                D = diag(1 ./ (Beta - sum((1 - opts.a2 * (Y .^ 2)) .* ex)));
                B = B + myy * B * (Y' * gauss - diag(Beta)) * D;

            % skew
            case {40}
                B = (X * ((X' * B) .^ 2)) / numwsmp;
            case {41}
                % Optimoitu
                Y = X' * B;
                Gskew = Y .^ 2;
                Beta = sum(Y .* Gskew);
                D = diag(1 ./ (Beta));
                B = B + myy * B * (Y' * Gskew - diag(Beta)) * D;
            case {42}
                Xsub = X(:, fi_getsmp(numwsmp, opts.smpsize));
                B = (Xsub * ((Xsub' * B) .^ 2)) / size(Xsub, 2);
            case {43}
                % Uusi optimoitu
                Y = X(:, fi_getsmp(numwsmp, opts.smpsize))' * B;
                Gskew = Y .^ 2;
                Beta = sum(Y .* Gskew);
                D = diag(1 ./ (Beta));
                B = B + myy * B * (Y' * Gskew - diag(Beta)) * D;
        end
    end

    if ~solutionfound
        warning( ...
            'neuroelf:FastICAWarn', ...
            'No convergence after %d steps.', ...
            maxiter ...
        );
        if ~isempty(B)
            % Symmetric orthogonalization.
            B = B * real(inv(B' * B)^(1/2));
            W = B' * whiteningMatrix;
            A = dewhiteningMatrix * B;
        else
            W = [];
            A = [];
        end
        if nargout == 2
            varargout{1} = A;
            varargout{2} = W;
        else
            if ~isempty(W)
                icasig = W * mixedsig + (W * mixedmean) * ones(1, nsmp);
            else
                icasig = [];
            end
            varargout{1} = icasig;
            if nargout > 1
                varargout{2} = A;
                varargout{3} = W;
            end
        end
        return;
    end

    % calculate ICA filters
    W = B' * whiteningMatrix;

% deflation
elseif opts.approach == 2

    % create array
    B = zeros(numvecs);

    % the search for a basis vector is repeated numics times
    round = 1;
    numFailures = 0;

    % loop for number of desired ICs
    while round <= opts.numics

        % progress
        if ~isempty(pbar)
            pbar.Progress(pstep * round, sprintf('Finding component %d...', round));
        end
        myy = myyOrig;
        usedNlinearity = gOrig;
        stroke = 0;
        notFine = true;
        long = false;
        endFinetuning = 0;

        % take a random initial vector of lenght 1 and orthogonalize it
        % with respect to the other vectors
        if ~opts.istate
            w = randn(numvecs, 1);
        else
            w = whiteningMatrix * opts.iguess(:, round);
        end

        % some computations
        w = w - B * B' * w;
        w = w / norm(w);

        wOld = zeros(size(w));
        wOld2 = zeros(size(w));

        % this is the actual fixed-point iteration loop
        iter = 1;
        gabba = 1;
        while iter <= opts.maxiter + gabba

            % project the vector into the space orthogonal to the space
            % spanned by the earlier found basis vectors; note that we can do
            % the projection with matrix B, since the zero entries do not
            % contribute to the projection
            w = w - B * B' * w;
            w = w / norm(w);

            if notFine
                if iter == opts.maxiter + 1
                    round = round - 1;
                    numFailures = numFailures + 1;
                    if numFailures > failureLimit
                        error( ...
                            'neuroelf:FastICAError', ...
                            'Unrecoverable error finding ICs.' ...
                        );
                    end
                    break;
                end
            else
                if iter >= endFinetuning
                    % so the algorithm will stop on the next test...
                    wOld = w;
                end
            end

            % test for termination condition; note that the algorithm has
            % converged if the direction of w and wOld is the same, this
            % is why we test the two cases
            if norm(w - wOld) < opts.epsilon || ...
                norm(w + wOld) < opts.epsilon
                if finetune && ...
                    notFine
                    notFine = false;
                    gabba = opts.maxftune;
                    wOld = zeros(size(w));
                    usedNlinearity = gFine;
                    myy = myyK * myyOrig;
                    endFinetuning = opts.maxftune + iter;

                else
                    numFailures = 0;
                    % Save the vector
                    B(:, round) = w;

                    % Calculate the de-whitened vector.
                    A(:, round) = dewhiteningMatrix * w;
                    % Calculate ICA filter.
                    W(round, :) = w' * whiteningMatrix;

                    % IC ready - next...
                    break;
                end

            elseif opts.stab
                if stroke == 0 && ...
                   (norm(w - wOld2) < opts.epsilon || ...
                    norm(w + wOld2) < opts.epsilon)
                    stroke = myy;
                    myy = 0.5 * myy;
                    if mod(usedNlinearity, 2) == 0
                        usedNlinearity = usedNlinearity + 1;
                    end
                elseif stroke > 0
                    myy = stroke;
                    stroke = 0;
                    if myy == 1 && ...
                        mod(usedNlinearity, 2) ~= 0
                        usedNlinearity = usedNlinearity - 1;
                    end
                elseif notFine && ...
                   ~long && ...
                    iter > (0.5 * opts.maxiter)
                    long = 1;
                    myy = 0.5 * myy;
                    if mod(usedNlinearity, 2) == 0
                        usedNlinearity = usedNlinearity + 1;
                    end
                end
            end

            % history of estimates
            wOld2 = wOld;
            wOld = w;

            % which sub-approach
            switch usedNlinearity

                % pow3
                case {10}
                    w = (X * ((X' * w) .^ 3)) / numwsmp - 3 * w;
                case {11}
                    EXGpow3 = (X * ((X' * w) .^ 3)) / numwsmp;
                    Beta = w' * EXGpow3;
                    w = w - myy * (EXGpow3 - Beta * w) / (3 - Beta);
                case {12}
                    Xsub = X(:, fi_getsmp(numwsmp, opts.smpsize));
                    w = (Xsub * ((Xsub' * w) .^ 3)) / size(Xsub, 2) - 3 * w;
                case {13}
                    Xsub = X(:, fi_getsmp(numwsmp, opts.smpsize));
                    EXGpow3 = (Xsub * ((Xsub' * w) .^ 3)) / size(Xsub, 2);
                    Beta = w' * EXGpow3;
                    w = w - myy * (EXGpow3 - Beta * w) / (3 - Beta);

                % tanh
                case {20}
                    if opts.a1 == 1
                        hypTan = fi_tanh(X' * w);
                    else
                        hypTan = fi_tanh(opts.a1 * (X' * w));
                    end
                    w = (X * hypTan - opts.a1 * sum(1 - hypTan .^ 2)' * w) / numwsmp;
                case {21}
                    if opts.a1 == 1
                        hypTan = fi_tanh(X' * w);
                    else
                        hypTan = fi_tanh(opts.a1 * (X' * w));
                    end
                    Beta = w' * X * hypTan;
                    w = w - myy * ((X * hypTan - Beta * w) / ...
                        (opts.a1 * sum((1 - hypTan .^ 2)') - Beta));
                case {22}
                    Xsub = X(:, fi_getsmp(numwsmp, opts.smpsize));
                    if opts.a1 == 1
                        hypTan = fi_tanh(Xsub' * w);
                    else
                        hypTan = fi_tanh(opts.a1 * (Xsub' * w));
                    end
                    w = (Xsub * hypTan - opts.a1 * ...
                        sum(1 - hypTan .^ 2)' * w) / size(Xsub, 2);
                case {23}
                    Xsub = X(:, fi_getsmp(numwsmp, opts.smpsize));
                    if opts.a1 == 1
                        hypTan = fi_tanh(Xsub' * w);
                    else
                        hypTan = fi_tanh(opts.a1 * (Xsub' * w));
                    end
                    Beta = w' * Xsub * hypTan;
                    w = w - myy * ((Xsub * hypTan - Beta * w) / ...
                        (opts.a1 * sum((1 - hypTan .^ 2)') - Beta));

                % gauss
                case {30}
                    % this has been split for performance reasons
                    u = X' * w;
                    u2 = u .^ 2;
                    ex = exp(-opts.a2 * u2 / 2);
                    gauss =  u .* ex;
                    dGauss = (1 - opts.a2 * u2) .* ex;
                    w = (X * gauss - sum(dGauss)' * w) / numwsmp;
                case {31}
                    u = X' * w;
                    u2 = u .^ 2;
                    ex = exp(-opts.a2 * u2 / 2);
                    gauss =  u .* ex;
                    dGauss = (1 - opts.a2 * u2) .* ex;
                    Beta = w' * X * gauss;
                    w = w - myy * ((X * gauss - Beta * w) / ...
                        (sum(dGauss)' - Beta));
                case {32}
                    Xsub = X(:, fi_getsmp(numwsmp, opts.smpsize));
                    u = Xsub' * w;
                    u2 = u .^ 2;
                    ex = exp(-opts.a2 * u2 / 2);
                    gauss =  u .* ex;
                    dGauss = (1 - opts.a2 * u2) .* ex;
                    w = (Xsub * gauss - sum(dGauss)' * w) / size(Xsub, 2);
                case {33}
                    Xsub = X(:, fi_getsmp(numwsmp, opts.smpsize));
                    u = Xsub' * w;
                    u2 = u .^ 2;
                    ex = exp(-opts.a2 * u2 / 2);
                    gauss =  u .* ex;
                    dGauss = (1 - opts.a2 * u2) .* ex;
                    Beta = w' * Xsub * gauss;
                    w = w - myy * ((Xsub * gauss - Beta * w) / ...
                        (sum(dGauss)' - Beta));
                % skew
                case {40}
                    w = (X * ((X' * w) .^ 2)) / numwsmp;
                case {41}
                    EXGskew = (X * ((X' * w) .^ 2)) / numwsmp;
                    Beta = w' * EXGskew;
                    w = w - myy * (EXGskew - Beta * w) / (-Beta);
                case {42}
                    Xsub = X(:, fi_getsmp(numwsmp, opts.smpsize));
                    w = (Xsub * ((Xsub' * w) .^ 2)) / size(Xsub, 2);
                case {43}
                    Xsub = X(:, fi_getsmp(numwsmp, opts.smpsize));
                    EXGskew = (Xsub * ((Xsub' * w) .^ 2)) / size(Xsub, 2);
                    Beta = w' * EXGskew;
                    w = w - myy * (EXGskew - Beta * w) / (-Beta);
            end

            % Normalize the new w.
            w = w / norm(w);
            iter = iter + 1;
        end
        round = round + 1;
    end
end

% in the end let's check the data for some security
if ~isreal(A)
  A = real(A);
  W = real(W);
end

% check for valid return
if ~isempty(W)
    icasig = W * mixedsig + (W * mixedmean) * ones(1, nsmp);
else
    icasig = [];
end

% the output depends on the number of output parameters
if nargout == 2
    varargout{1} = A;
    varargout{2} = W;
else
    varargout{1} = icasig;
    if nargout > 1
        varargout{2} = A;
        varargout{3} = W;
    end
end

% delete progress bar
if ~isempty(pbar) && ...
    isempty(opts.pbar)
    closebar(pbar);
end



% subfunctions



% calculate tanh simplier and faster than Matlab tanh
function y = fi_tanh(x)
y = 1 - 2 ./ (exp(2 * x) + 1);

% create (boolean) index vector to retrieve a smaller number of samples
function smp = fi_getsmp(max, percentage)
smp = (rand(1, max) < percentage);

% calculate PCA for data
function [E, D] = fi_pcamat(vectors, eig1, eign)

% check the optional parameters
oldDimension = size(vectors, 1);

% calculate the covariance matrix
covarianceMatrix = cov(vectors', 1);

% calculate the eigenvalues and eigenvectors of covariance matrix
[E, D] = eig(covarianceMatrix);

% the rank is determined from the eigenvalues - and not directly by
% using the function rank - because function rank uses svd, which
% in some cases gives a higher dimensionality than what can be used
% with eig later on (eig then gives negative eigenvalues).
rankTolerance = 1e-7;
maxeign = sum(diag(D) > rankTolerance);
if maxeign == 0,
    error( ...
        'neuroelf:FastICAError', ...
        'All eigenvalues are smaller than the tolerance.' ...
    );
end

% sort the eigenvalues - decending
eigenvalues = sort(diag(D), 'descend');

% restrict to actual maximum
eig1 = min(eig1, maxeign);
eign = min(eign, maxeign);

% drop the smaller eigenvalues
if eign < oldDimension
    lowerLimitValue = (eigenvalues(eign) + eigenvalues(eign + 1)) / 2;
else
    lowerLimitValue = eigenvalues(oldDimension) - 1;
end
lowerColumns = diag(D) > lowerLimitValue;

% drop the larger eigenvalues
if eig1 > 1
    higherLimitValue = (eigenvalues(eig1 - 1) + eigenvalues(eig1)) / 2;
else
    higherLimitValue = eigenvalues(1) + 1;
end
higherColumns = diag(D) < higherLimitValue;

% combine the results from above
selectedColumns = lowerColumns & higherColumns;

% some last check
if sum(selectedColumns) ~= (eign - eig1 + 1)
    error( ...
        'neuroelf:FastICAError', ...
        'Selected a wrong number of dimensions.' ...
    );
end

% select the colums which correspond to the desired range of eigenvalues
E = E(:, selectedColumns);
D = D(selectedColumns, selectedColumns);


% fi_remmean - remove the mean from vectors
function [newVectors, meanValue] = fi_remmean(vectors)
meanValue = mean(vectors, 2);
newVectors = vectors - meanValue * ones (1, size(vectors, 2));


% fi_whitenv - whitenv vectors
function [newVectors, whiteningMatrix, dewhiteningMatrix] = fi_whitenv(vectors, E, D)
% in some cases, rounding errors in Matlab cause negative eigenvalues
% (elements in the diagonal of D); since it is difficult to know when
% this happens, it is difficult to correct it automatically-- instead
% an error is raised (for now...)
if any (diag (D) < 0)
    error( ...
        'neuroelf:FastICAError', ...
        '%d negative eigenvalues computed from the cov matrix.', ...
        sum(diag(D) < 0) ...
    );
end

% calculate the whitening and dewhitening matrices (these handle
% dimensionality simultaneously)
whiteningMatrix = inv(sqrt(D)) * E';
dewhiteningMatrix = E * sqrt(D);

% project to the eigenvectors of the covariance matrix
% whiten the samples and reduce dimension simultaneously.
newVectors =  whiteningMatrix * vectors;

% as security measure...
if ~isreal(newVectors)
    error( ...
        'neuroelf:FastICAError', ...
        'Whitened vectors have imaginary values.' ...
    );
end
