hqusers_num = 150;
lqusers_num = 150;
mausers_num_max = 100;

tasks_num = 1600;

hqusers_mat = zeros(hqusers_num, tasks_num);
for i=1:hqusers_num
    ah = rand()*5 + 10; % random number in the range [10-15]
    bh = rand()*2 + 3; % random number in the range [3.0-5.0]
    hqusers_pd = makedist('Beta', 'a', ah, 'b', bh);
    hqusers_mat(i,:) = random(hqusers_pd, 1, tasks_num);
end

lqusers_mat = zeros(lqusers_num, tasks_num);
for i=1:lqusers_num
    al = rand()*3 + 9; % random number in the range [9-12]
    bl = rand()*2 + 7; % random number in the range [7-9]
    lqusers_pd = makedist('Beta', 'a', al, 'b', bl);
    lqusers_mat(i,:) = random(lqusers_pd, 1, tasks_num);   
end

mausers1_mat = zeros(mausers_num_max, tasks_num);
mausers2_mat = zeros(mausers_num_max, tasks_num);
mausers_mat_max = zeros(mausers_num_max, tasks_num);
for i=1:mausers_num_max
    % High quality beta distribution model
    am1 = rand()*4 + 18; % random number in the range [18-22]
    bm1 = rand()*1 + 2.5; % random number in the range [2.5-3.5]
    amusers_pd1 = makedist('Beta', 'a', am1, 'b', bm1);
    mausers1_mat(i,:) = random(amusers_pd1, 1, tasks_num);
    % Low quality beta distribution model
    am2 = rand()*2 + 4; % random number in the range [4-6]
    bm2 = rand()*10 + 25; % random number in the range [25-35]
    amusers_pd2 = makedist('Beta', 'a', am2, 'b', bm2);
    mausers2_mat(i,:) = random(amusers_pd2, 1, tasks_num);    

    % the PDF for the bimodal distribution of the two beta distribution with the alpha mixture coefficient
    % bimodal_betapdf = (1-alpha)*y3 + alpha*y4
    alpha = 0.7; % mixture coefficient
    for j=1:tasks_num
        r = rand();
        if (r < alpha)
            mausers_mat_max(i,j) = mausers1_mat(i,j);
        else
            mausers_mat_max(i,j) = mausers2_mat(i,j);
        end
    end
end


provider_min = 5;
%provider_max = round(users_num/2);
iter_max = 100;
result = zeros(iter_max, 4);

for iter=1:iter_max
    mausers_num = iter;
    mausers_mat = mausers_mat_max(1:iter, :);
    
    users_num = mausers_num + hqusers_num + lqusers_num;
    provider_max = round(users_num/2);    
    
    % Generate a trust relationship matrix between 100 users
    % Trust scores at the bootstrap are equally set as 0.5
    trust_init_value = 0.5;
    exp_init_value = 0.5;
    qod_init_value = 0.5;
    trust_mat = ones(users_num)*trust_init_value;
    exp_mat = ones(users_num)*exp_init_value;

    interaction_trust_mat = zeros(users_num);
    interaction_qod_mat = zeros(users_num);
    interaction_avg_mat = zeros(users_num);
    interaction_reg_mat = zeros(users_num);

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
    trust_based_qod = zeros(tasks_num, 1); % Average QoD for the Highest Trust-based MCS scheme (Trust-based scheme)
    qod_based_qod = zeros(tasks_num, 1); % Average QoD for the Highest Previous QoD scheme (QoD-based Scheme)
    avg_based_qod = zeros(tasks_num, 1); % Average QoD for the Highest Average QoD scheme (Average-based Scheme)
    qos_reg_based = zeros(tasks_num, 1); % QoS record for the Highest Average QoD scheme (Regression Scheme)

    % Create a matrix of QoD from all users in all tasks.
    A = [hqusers_mat; lqusers_mat; mausers_mat];
    % Random Permutation for Matrix A
    shuffledA = A(randperm(size(A,1)),:);

    % For a QoD-based scheme, the initial step is added in which all users have the same QoD at 0.5
    A1 = [ones(users_num, 1)*qod_init_value shuffledA]; % For the initial step: no previous task, QoD of all users are same at 0.5
    % For Average-based scheme, initial step sets all users a same average of QoD at 0.5
    Avg = ones(users_num, 1)*qod_init_value;

     % For the Regression scheme, the initial step is added in which all users have the same QoD
    QoD_mat_reg_scheme = [ones(users_num, 4).*[0.5 0.55 0.65 0.6] shuffledA]; % For the initial step: no previous task, QoD of all users are same at 0.8

    for i=1:tasks_num
        % Randomly select a user as the service requester
        requester = round(rand()*(users_num-1) + 1); % random number in the range [1-users_num]
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

        % Calculate the QoS for the task i
        trust_based_qod(i) = provider_num/abs(log(prod(shuffledA(maxValueIndices_trust, i)))); % take the average

        
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
        % Calculate the QoS for the task i
        qod_based_qod(i) = provider_num/abs(log(prod(A1(maxValueIndices_qod, i+1))));
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Average-based MCS System %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Avg_temp = Avg;
        Avg_temp(requester) = 0;
        [sortedX_avg,sortingIndices_avg] = sort(Avg_temp,'descend');
        maxValueIndices_avg = sortingIndices_avg(1:provider_num);

        % Update interaction matrix
        interaction_avg_mat(requester, maxValueIndices_avg) = interaction_qod_mat(requester, maxValueIndices_avg) + 1; 
        % Calculate the QoS for the task i
        avg_based_qod(i) = provider_num/abs(log(prod(shuffledA(maxValueIndices_avg, i))));
        
        % Update Avg
        Avg(maxValueIndices_avg) = ((Avg(maxValueIndices_avg) .* (sum(interaction_avg_mat(:, maxValueIndices_avg)).') + shuffledA(maxValueIndices_avg, i)))...
            ./((sum(interaction_avg_mat(:, maxValueIndices_avg))).' + 1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Regression MCS System %
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        Reg_QoD_mat = QoD_mat_reg_scheme(:, 1:i+3); % getting data for regression model
        Reg_QoD_mat(requester, :) = 0; % requester data should be clear in order not to be chosen
        t = linspace(1,i+3,i+3);
        f = zeros(users_num, 1);
        for k=1:users_num
            fn = polyfit(t, Reg_QoD_mat(k,:), 3); % figure out the coefficients
            f(k) = polyval(fn, i+4); % predict the next QoD scores
        end
        [sortedX_reg,sortingIndices_reg] = sort(f,'descend'); % find the users with highest predicted QoD scores
        maxValueIndices_reg = sortingIndices_reg(1:provider_num);
                
        % Update interaction matrix
        interaction_reg_mat(requester, maxValueIndices_reg) = interaction_reg_mat(requester, maxValueIndices_reg) + 1; 
        % Calculate the QoS for the task i
        qos_reg_based(i) = provider_num/abs(log(prod(QoD_mat_reg_scheme(maxValueIndices_reg, i+4))));
    end
    result(iter, 1) = sum(trust_based_qod)/tasks_num;
    result(iter, 2) = sum(qod_based_qod)/tasks_num;
    result(iter, 3) = sum(avg_based_qod)/tasks_num;
    result(iter, 4) = sum(qos_reg_based)/tasks_num;
end
save result_qos_fig2_1600tasks_mauser_varied_3degree
