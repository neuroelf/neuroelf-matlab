function [plab, acc, dvpe] = ne_svmpredict(labels, tdata, model, opts)
% ne_svmpredict  - predict labels using libSVM based classifier
%
% FORMAT:       [plab, acc, dvpe] = ne_svmpredict(labels, tdata, model, opts)
%
% Input fields:
%
%       labels      Lx1 double array with labels
%       tdata       CxT double testing case data (with T features/dims)
%       model       1x1 struct (output of ne_svmtrain)
%       opts        libSVM options
%        .probest   flag, prob. estimates for SVC/SVR (default: false)
%        .quiet     libSVM quiet mode (default: true)
%        .version   libSVM version to use, either 2 or {3}
%
% Output fields:
%
%       plab        predicted labels for testing case data
%       acc         3x1 accuracy, where
%           (1)     percent accurate (0...100)
%           (2)
%           (3)
%       dvpe        decision values/probability estimates
%
% Note: this functions uses the MEX-file svmtrainc, which is courtesy
%       of Chih-Jen Lin, being part of the the libSVM library,
%       see http://www.csie.ntu.edu.tw/~cjlin/libsvm/ for details

% Version:  v1.1
% Build:    16032111
% Date:     Mar-21 2016, 11:34 AM EST
% Author:   Chih-Jen-Lin, Department of CS, National Taiwan University
% Editor:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://www.csie.ntu.edu.tw/~cjlin/libsvm/

% rudimentary argument check
if nargin < 3 || ...
   ~isa(labels, 'double') || ...
    isempty(labels) || ...
    numel(labels) ~= size(labels, 1) || ...
    any(isinf(labels) | isnan(labels)) || ...
   ~isa(tdata, 'double') || ...
    isempty(tdata) || ...
    size(tdata, 1) ~= numel(labels) || ...
    ndims(tdata) > 2 || ...
    any(isinf(tdata(:)) | isnan(tdata(:))) || ...
   ~isstruct(model) || ...
    numel(model) ~= 1 || ...
   ~isfield(model, 'Parameters') || ...
   ~isfield(model, 'Label') || ...
   ~isfield(model, 'sv_coef') || ...
   ~isfield(model, 'SVs')
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing argument in call to ne_svmpredict.' ...
    );
end

% allow libSVM string directly
if nargin > 3 && ...
    ischar(opts) && ...
    numel(opts) == size(opts, 2)
    optsstring = opts;

    % make version 3 default
    opts = struct('version', 3);

% otherwise parse/create options struct
else

    % create struct if necessary
    if nargin < 4 || ...
       ~isstruct(opts) || ...
        numel(opts) ~= 1
        opts = struct;
    end

    % parse individual options
    if ~isfield(opts, 'probest') || ...
       ~islogical(opts.probest) || ...
        numel(opts.probest) ~= 1 || ...
       ~opts.probest
        opts.probest = 0;
    else
        opts.probest = 1;
    end
    if ~isfield(opts, 'quiet') || ...
       ~islogical(opts.quiet) || ...
        numel(opts.quiet) ~= 1 || ...
        opts.quiet
        quietstr = '-q ';
    else
        quietstr = '';
    end
    if ~isfield(opts, 'version') || ...
       ~isa(opts.version, 'double') || ...
        numel(opts.version) ~= 1 || ...
        isinf(opts.version) || ...
        isnan(opts.version) || ...
       ~any(opts.version == [2, 3])
        opts.version = 3;
    end

    % combine to libsvm-compatible string
    optsstring = sprintf('%s-b %d', quietstr, opts.probest);
end

% pass on
if opts.version < 3
    [plab, acc, dvpe] = svmpredictc2x(labels, tdata, model, optsstring);
else
    [plab, acc, dvpe] = svmpredictc3x(labels, tdata, model, optsstring);
end
