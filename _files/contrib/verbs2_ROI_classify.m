% number of classify cross-validation folds
nfolds = 5000;

% number of hold-out samples (per class)
nholdout = 1;

% use NeuroElf
n = neuroelf;

% load GLM
if exist('glm', 'var') == 0 || ...
   ~isxff(glm, 'glm')
    glm = bless(xff('E:/verbs2/Analysis_SPM5/GLM_st-psg_FFX.glm'));
end

% load VOI
if exist('voi', 'var') == 0 || ...
   ~isxff(voi, 'voi')
    voi = bless(xff('*.voi', 'Please select VOI to use as classification region...'));
end

% access (individual voxel) raw data
[vb, vbv] = glm.VOIBetas(voi, struct('vl', 1));

% number of subjects
subs = glm.Subjects;
nsubs = numel(subs);

% get subject predictors
sps = glm.SubjectPredictors;

% get SensNeg and ResNeg entries
sensneg = find(~cellfun('isempty', regexpi(sps, '^Pic_Sens')));
resneg = find(~cellfun('isempty', regexpi(sps, '^Pic_Res')));

% total number of elements
e1 = numel(sensneg);
e2 = numel(resneg);
f1 = e1 - nholdout;
f2 = e2 - nholdout;
t1 = f1 + 1;
t2 = f2 + 1;

% training labels
tl = [ones(f1, 1); 2 .* ones(f2, 1)];
tsl = [ones(nholdout, 1); 2 .* ones(nholdout, 1)];

% generate matrix holding the N-fold accuracy values
class_accuracy = NaN .* zeros(nsubs, nfolds);

% iterate over subjects
for sc = 1:nsubs
    
    % show progress
    fprintf('Working on subject %s...\n', subs{sc});
    pause(0.01);
    
    % create random vectors
    [rvd, rvi1] = sort(randn(e1, nfolds));
    [rvd, rvi2] = sort(randn(e2, nfolds));
    
    % iterate over folds
    for fc = 1:nfolds
        
        % concatenate data over selected folds
        d1 = double(cat(2, vbv{sc, sensneg(rvi1(1:f1, fc))})');
        d2 = double(cat(2, vbv{sc, resneg(rvi2(1:f1, fc))})');
        
        % classifier training
        trained = n.ne_svmtrain(tl, [d1; d2]);
        
        % concatenate test data
        testdata = cat(1, ...
            double(cat(2, vbv{sc, sensneg(rvi1(t1:end, fc))})'), ...
            double(cat(2, vbv{sc, resneg(rvi2(t2:end, fc))})'));
        
        % check accuracy
        [plab, acc] = n.ne_svmpredict(tsl, testdata, trained);
        
        % store accuracy
        class_accuracy(sc, fc) = 0.01 * acc(1);
    end
end
