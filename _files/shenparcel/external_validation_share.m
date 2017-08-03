% Copyright 2015 Monica Rosenberg, Emily Finn, and Dustin Scheinost

% This code is released under the terms of the GNU GPL v2. This code
% is not FDA approved for clinical use; it is provided
% freely for research purposes. If using this in a publication
% please reference this properly as: 

% Rosenberg MD, Finn ES, Scheinost D, Papademetris X, Shen X, 
% Constable RT & Chun MM. (2016). A neuromarker of sustained
% attention from whole-brain functional connectivity. Nature Neuroscience 
% 19(1), 165-171.

% This code provides a framework for implementing functional
% connectivity-based behavioral prediction across datasets, as 
% described in Rosenberg, Finn et al 2016 (see above for full reference). 
% The first input ('train_mats') is a pre-calculated MxMxN matrix 
% containing all individual-subject connectivity matrices in the training 
% set, where M = number of nodes in the chosen brain atlas and N = number of
% subjects. Each element (i,j,k) in these matrices represents the
% correlation between the BOLD timecourses of nodes i and j in subject k
% during a single fMRI session. The second input ('behav') is the
% Nx1 vector of scores for the behavior of interest for all subjects.
% The input 'validation_mats' is a pre-calculated MxMxP matrix 
% containing all individual-subject connectivity matrices in the test set, 
% where P = number of validaiton subjects in the test set. 'validation_behav' 
% is a Px1 vector of scores for the behavior of interest in the test set. 

% As in the reference paper, the predictive power of the model is assessed
% via correlation between predicted and observed scores across all
% subjects. Note that this assumes normal or near-normal distributions for
% both vectors, and does not assess absolute accuracy of predictions (only
% relative accuracy within the sample). It is recommended to explore
% additional/alternative metrics for assessing predictive power, such as
% prediction error sum of squares or prediction r^2.

clear;
clc;

overlap = 0;                                % overlap = 1 if you want to select edges for predictive networks if they appear in every round of leave-one-out cross-validation in the training set. 
                                            % overlap = 0 if you want to define predictive networks using all training subjects (i.e., no leave-one-out in the training set)
n_node  = 268;                              % number of nodes
thresh  = 0.01;                             % threshold for feature selection

% Training data 
train_mats  = ;                             % training data (n_node x n_node x n_sub symmetrical connectivity matrices)
behav       = ;                             % n_sub x 1 vector of behavior
n_sub       = size(train_mats,3);           % total number of training subjects
n_train_sub = n_sub-1;                      % number of training subjects in each round of leave-one-out

% Validation data
validation_mats  = ;                        % validation data (n_node x n_node x n_validation_sub symmetrical connectivity matrices)
validation_behav = ;                        % n_validation_sub x 1 vector of behavior
n_validation_sub = size(validation_mats,3); % total number of validation subjects

aa     = ones(n_node,n_node);
aa_upp = triu(aa,1);
upp_id = find(aa_upp);                      % indices of edges in the upper triangular of an n_node x n_node matrix
n_edge = length(upp_id);                    % total number of edges

if overlap == 1
    
    pos_mask_all = zeros(n_node, n_node, n_sub);
    neg_mask_all = zeros(n_node, n_node, n_sub);
    
    for excl_sub = 1:n_sub; % create network masks for every iteration of leave-one-subject-out cross-validation in the training set

        excl_sub
        
        % exclude data from left-out subject
        train_mats_tmp = train_mats;
        train_mats_tmp(:,:,excl_sub) = [];
        train_behav = behav;
        train_behav(excl_sub) = [];
        
        % create n_train_sub x n_edge matrix
        train_vect = reshape(train_mats_tmp, n_node*n_node, n_train_sub)';
        upp_vect   = train_vect(:,upp_id);
        
        % relate behavior to edge strength across training subjects
        cp = zeros(n_edge, 1);
        cr = zeros(n_edge, 1);
        
        for ii = 1:n_edge
            [b,stats] = robustfit(upp_vect(:,ii), train_behav);
            cp(ii)    = stats.p(2);
            cr(ii)    = sign(stats.t(2))*sqrt((stats.t(2)^2/(n_train_sub-2))/(1+(stats.t(2)^2/(n_train_sub-2))));
        end
        
        % select edges based on threshold
        pos_edge = zeros(1, n_edge);
        neg_edge = zeros(1, n_edge);
        
        cp_pos           = find(cp<thresh & cr>0);
        pos_edge(cp_pos) = 1;
        cp_neg           = find(cp<thresh & cr<0);
        neg_edge(cp_neg) = 1;
        
        pos_mask = zeros(n_node, n_node);
        neg_mask = zeros(n_node, n_node);
        
        pos_mask(upp_id) = pos_edge;  % Here, masks are NOT symmetrical. To make symmetrical, set pos_mask = pos_mask + pos_mask'
        neg_mask(upp_id) = neg_edge;
        
        % save masks from every iteration of the leave-one-out loop
        pos_mask_all(:,:,excl_sub) = pos_mask;
        neg_mask_all(:,:,excl_sub) = neg_mask;
    end
    
    % select edges that appear in all iteraitons of leave-one-out
    pos_overlap = zeros(n_node, n_node);
    neg_overlap = zeros(n_node, n_node);
    
    pos_overlap(sum(pos_mask_all,3) == n_sub) = 1;
    neg_overlap(sum(neg_mask_all,3) == n_sub) = 1;
    
elseif overlap == 0
        
    % create n_sub x n_edge matrix
    train_vect = reshape(train_mats, n_node*n_node, n_sub)';
    upp_vect   = train_vect(:,upp_id);

    % relate behavior to edge strength across ALL subjects in training set
    cp = zeros(n_edge, 1);
    cr = zeros(n_edge, 1);

    for ii = 1:n_edge
        [b,stats] = robustfit(upp_vect(:,ii), behav);
        cp(ii)    = stats.p(2);
        cr(ii)    = sign(stats.t(2))*sqrt((stats.t(2)^2/(n_sub-2))/(1+(stats.t(2)^2/(n_sub-2))));
    end

    % select edges based on threshold
    pos_edge = zeros(1, n_edge);
    neg_edge = zeros(1, n_edge);

    cp_pos           = find(cp<thresh & cr>0);
    pos_edge(cp_pos) = 1;
    cp_neg           = find(cp<thresh & cr<0);
    neg_edge(cp_neg) = 1;

    pos_mask = zeros(n_node, n_node);
    neg_mask = zeros(n_node, n_node);

    pos_mask(upp_id) = pos_edge; % Here, masks are NOT symmetrical. To make symmetrical, set pos_mask = pos_mask + pos_mask'
    neg_mask(upp_id) = neg_edge;

    pos_overlap = pos_mask;
    neg_overlap = neg_mask;
end

% sum edges for all subjects in the training set
train_pos_sum = zeros(n_sub,1);
train_neg_sum = zeros(n_sub,1);
    
for k = 1:n_sub
    train_pos_sum(k) = sum(sum(pos_overlap.*train_mats(:,:,k)));
    train_neg_sum(k) = sum(sum(neg_overlap.*train_mats(:,:,k)));
end
    
% build model with training data
b_pos      = robustfit(train_pos_sum, behav);
b_neg      = robustfit(train_neg_sum, behav);
robGLM_fit = robustfit([train_pos_sum train_neg_sum],behav);
    
% generate predictions for validation set
pred_pos = zeros(n_validation_sub,1);
pred_neg = zeros(n_validation_sub,1);
pred_glm = zeros(n_validation_sub,1);

validation_pos_sum = zeros(n_validation_sub,1);
validation_neg_sum = zeros(n_validation_sub,1);

for vs = 1:n_validation_sub
    validation_pos_sum(vs) = sum(sum(pos_overlap.*validation_mats(:,:,vs)));
    validation_neg_sum(vs) = sum(sum(neg_overlap.*validation_mats(:,:,vs)));

    pred_pos(vs) = (b_pos(2)*validation_pos_sum(vs)) + b_pos(1);
    pred_neg(vs) = (b_neg(2)*validation_neg_sum(vs)) + b_neg(1);
    pred_glm(vs) = robGLM_fit(1) + robGLM_fit(2)*validation_pos_sum(vs) + robGLM_fit(3)*validation_neg_sum(vs);
end

% correlate predicted and observed behavior
[r_pos, p_pos] = corr(validation_behav, pred_pos);
[r_neg, p_neg] = corr(validation_behav, pred_neg);
[r_glm, p_glm] = corr(validation_behav, pred_glm);

ax1 = subplot(2,2,1);
scatter(validation_behav, pred_pos)
title(['pos r = ' num2str(round(r_pos*100)/100) ', p = ' num2str(p_pos)])
xlabel('Observed')
ylabel('Predicted')

ax2 = subplot(2,2,2);
scatter(validation_behav, pred_neg)
title(['neg r = ' num2str(round(r_neg*100)/100) ', p = ' num2str(p_neg)])
xlabel('Observed')
ylabel('Predicted')

ax3 = subplot(2,2,3);
scatter(validation_behav, pred_glm)
title(['glm r = ' num2str(round(r_glm*100)/100) ', p = ' num2str(p_glm)])
xlabel('Observed')
ylabel('Predicted')