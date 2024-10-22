hqusers_num = 100;
lqusers_num = 150;
mausers_num = 50;
users_num = mausers_num + hqusers_num + lqusers_num;
provider_max = round(users_num/2);

provider_min = 5;
iter_max = 101;
result = zeros(iter_max, 3);
tasks_num_max = 1600;

hqusers_mat_max = zeros(hqusers_num, tasks_num_max);
for i=1:hqusers_num
    ah = rand()*5 + 10; % random number in the range [10-15]
    bh = rand()*2 + 3; % random number in the range [3.0-5.0]
    hqusers_pd = makedist('Beta', 'a', ah, 'b', bh);
    hqusers_mat_max(i,:) = random(hqusers_pd, 1, tasks_num_max);
end

lqusers_mat_max = zeros(lqusers_num, tasks_num_max);
for i=1:lqusers_num
    al = rand()*3 + 9; % random number in the range [9-12]
    bl = rand()*2 + 7; % random number in the range [7-9]
    lqusers_pd = makedist('Beta', 'a', al, 'b', bl);
    lqusers_mat_max(i,:) = random(lqusers_pd, 1, tasks_num_max);   
end

mausers1_mat = zeros(mausers_num, tasks_num_max);
mausers2_mat = zeros(mausers_num, tasks_num_max);
mausers_mat_max = zeros(mausers_num, tasks_num_max);
for i=1:mausers_num
    % High quality beta distribution model
    am1 = rand()*4 + 18; % random number in the range [18-22]
    bm1 = rand()*1 + 2.5; % random number in the range [2.5-3.5]
    amusers_pd1 = makedist('Beta', 'a', am1, 'b', bm1);
    mausers1_mat(i,:) = random(amusers_pd1, 1, tasks_num_max);
    % Low quality beta distribution model
    am2 = rand()*2 + 4; % random number in the range [4-6]
    bm2 = rand()*10 + 25; % random number in the range [25-35]
    amusers_pd2 = makedist('Beta', 'a', am2, 'b', bm2);
    mausers2_mat(i,:) = random(amusers_pd2, 1, tasks_num_max);    

    % the PDF for the bimodal distribution of the two beta distribution with the alpha mixture coefficient
    % bimodal_betapdf = (1-alpha)*y3 + alpha*y4
    alpha = 0.7; % mixture coefficient
    for j=1:tasks_num_max
        r = rand();
        if (r < alpha)
            mausers_mat_max(i,j) = mausers1_mat(i,j);
        else
            mausers_mat_max(i,j) = mausers2_mat(i,j);
        end
    end
end


% Randomly select a requester from an array
% Create a list of requesters for a number of service requests
requester_max = round(rand(tasks_num_max, 1)*(users_num-1) + 1); % random number in the range [1-users_num]

for iter=1:iter_max    
    tasks_num = 100 + 15*(iter-1);    
    hqusers_mat = hqusers_mat_max(:, 1:tasks_num); 
    lqusers_mat = lqusers_mat_max(:, 1:tasks_num);
    mausers_mat = mausers_mat_max(:, 1:tasks_num);
    
    % Generate a trust relationship matrix between users
    % Trust scores at the bootstrap are equally set as 0.5
    trust_init_value = 0.5;
    exp_init_value = 0.5;
    qod_init_value = 0.5;
    trust_mat = ones(users_num)*trust_init_value;
    exp_mat = ones(users_num)*exp_init_value;

    interaction_trust_mat = zeros(users_num);
    interaction_qod_mat = zeros(users_num);
    interaction_avg_mat = zeros(users_num);

    % Set the diagonal elements = 0
    trust_mat(logical(eye(size(trust_mat)))) = 0;
    exp_mat(logical(eye(size(exp_mat)))) = 0;

    % For calculating Experience, it is needed to keep track of previous experience values+
    exp_prev_mat = exp_mat;
    % For calculating Reputation, divide into negative experiences and positive experience matrices
    rep_mat = ones(users_num, 1)/users_num;
    rep_pos = ones(users_num, 1)/users_num; % positive reputation
    rep_neg = ones(users_num, 1)/users_num; % negative reputation

    % Keep track of the average QoD for the MCS system for each task. There are tasks_num tasks so that we store the average QoD in the arrays:
    trust_based_qod = zeros(tasks_num, 1); % Average System Damage for the Highest Trust-based MCS scheme (Trust-based scheme)
    qod_based_qod = zeros(tasks_num, 1); % Average System Damage for the Highest Previous QoD scheme (QoD-based Scheme)
    avg_based_qod = zeros(tasks_num, 1); % Average System Damage for the Highest Average QoD scheme (Average-based Scheme)

    % Create a matrix of QoD from all users in all tasks.
    A = [hqusers_mat; lqusers_mat; mausers_mat];
    % Random Permutation for Matrix A
    shuffledA = A(randperm(size(A,1)),:);

    % For a QoD-based scheme, the initial step is added in which all users have the same QoD at 0.8
    A1 = [ones(users_num, 1)*qod_init_value shuffledA]; % For the initial step: no previous task, QoD of all users are same at 0.8
    % For Average-based scheme, initial step sets all users a same average of QoD at 0.5
    Avg = ones(users_num, 1)*qod_init_value;

    for i=1:tasks_num
        % Randomly select a user as the service requester
        requester = requester_max(i); % random number in the range [1-users_num]
        % Randomly select a number of users among the remaining users as the service providers 
        provider_num = round(rand()*(provider_max - provider_min) + provider_min); % random number in the range [provider_min, provider_max]

        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Trust-based MCS System %
        %%%%%%%%%%%%%%%%%%%%%%%%%%

        % Select the recruited users which have the highest trust relationship with the requester.    
        [sortedX,sortingIndices_trust] = sort(trust_mat(requester, :),'descend');
        % maxValues = sortedX(1:provider_num);
        maxValueIndices_trust = sortingIndices_trust(1:provider_num);

        for j=1:provider_num % Consider all QoD from provider_num contributors (data providers)        
            % Update Experience value between the requester and the selected user (maxValueIndices_trust(j))
            % function [ newExp ] = expCalculation( fbValue, currentExp, previousExp)
            newExp = expCalculation(shuffledA(maxValueIndices_trust(j), i), exp_mat(requester, maxValueIndices_trust(j)), exp_prev_mat(requester, maxValueIndices_trust(j)));        
            exp_prev_mat(requester, maxValueIndices_trust(j)) = exp_mat(requester, maxValueIndices_trust(j)); % Update previous_exp_mat
            exp_mat(requester, maxValueIndices_trust(j)) = newExp; % Update the exp_mat
        end

        % Update interaction matrix
        interaction_trust_mat(requester, maxValueIndices_trust) = interaction_trust_mat(requester, maxValueIndices_trust) + 1;

        % Update Reputation value    
        [rep_pos, rep_neg] = repCalculation(rep_pos, rep_neg, exp_mat);
        rep1 = rep_pos - rep_neg;    
        rep = rep1*users_num; 
        % rep(rep < 0) = 0;   
        %rep = zeros(1, users_num);

        % Update Trust value
        trust_mat = (exp_mat + rep.')/2;

        % Calculate the average QoD for the task i
        trust_based_qod(i) = log(prod(shuffledA(maxValueIndices_trust, i))); % take the average

        %%%%%%%%%%%%%%%%%%%%%%%%
        % QoD-based MCS System %
        %%%%%%%%%%%%%%%%%%%%%%%%

        % QoD-based MCS
        % Choose the users which provide the best QoD values except the requester
        A2 = A1(:, i); 
        A2(requester) = 0;

        [sortedX_qod,sortingIndices_qod] = sort(A2,'descend');
        maxValueIndices_qod = sortingIndices_qod(1:provider_num);
        % Update interactions matrix    
        interaction_qod_mat(requester, maxValueIndices_qod) = interaction_qod_mat(requester, maxValueIndices_qod) + 1;    
        % Calculate the average QoD for the task i
        qod_based_qod(i) = log(prod(A1(maxValueIndices_qod, i+1))); % take the average

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Average-based MCS System %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Avg_temp = Avg;
        Avg_temp(requester) = 0;
        [sortedX_avg,sortingIndices_avg] = sort(Avg_temp,'descend');
        maxValueIndices_avg = sortingIndices_avg(1:provider_num);

        % Update interaction matrix
        interaction_avg_mat(requester, maxValueIndices_avg) = interaction_qod_mat(requester, maxValueIndices_avg) + 1; 
        % Calculate the average QoD for the task i
        avg_based_qod(i) = log(prod(shuffledA(maxValueIndices_avg, i)));
        % Update Avg
        Avg(maxValueIndices_avg) = ((Avg(maxValueIndices_avg) .* (sum(interaction_avg_mat(:, maxValueIndices_avg)).') + shuffledA(maxValueIndices_avg, i)))...
            ./((sum(interaction_avg_mat(:, maxValueIndices_avg))).' + 1);
    end
    result(iter, 1) = sum(trust_based_qod)/tasks_num;
    result(iter, 2) = sum(qod_based_qod)/tasks_num;
    result(iter, 3) = sum(avg_based_qod)/tasks_num;
end
