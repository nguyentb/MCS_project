% Experience Model - Graph Drawing for paperwork
% Experience Model specifies three trends: Development, Loss and Decay
% Experience is normalized in the range [0,1]
%
%--------------------------------------------------------------------
% Randomly create a series of interactions
n = 70; % number of interactions
coop_threshold = 0.6;
uncoop_threshold = 0.3;
%interaction = rand(n, 1);
%coop_interaction = interaction(interaction > coop_threshold);
%uncoop_interaction = interaction(interaction < uncoop_threshold);
coop_interaction = rand(n, 1)*(1 - coop_threshold) + coop_threshold;
uncoop_interaction = rand(n, 1)*uncoop_threshold;

exp_increase = zeros(n,0);
exp_init = 0.3;
exp_increase(1) = exp_init;

grid on;
grid minor;
hold on;

% Development Model
% x(i+1) = x(i) + interaction_value*alpha*(1-x(i));
% alpha is a parameter defining the maximum increase value
alpha = 0.15; % Maximum increase value
for i=1:n-1
    exp_increase(i+1)= exp_increase(i) + coop_interaction(i)*alpha*(1-exp_increase(i));
end
plot(0:n-1, exp_increase, 'b-.*');

% Decay Model: let define alpha = maximum increase of Experience
% exp(t+1) = exp(t) - decay
%let set the parameters
min_decay = 0.005;
decay_rate = 0.005;
min_exp = 0;

exp_decay = zeros(n, n);
exp_decay(:, 1) = exp_increase;
exp_decay(:, 2) = exp_increase;

for i=2:n
    for j=2:n-1
        exp_decay(i, j+1) = max(exp_init, exp_decay(i, j) - min_decay*(1 + decay_rate - exp_decay(i, j-1)));    
    end
end
decay_start = 5;
plot(decay_start:n, exp_decay(decay_start,2:(n-decay_start + 2)), '-.', 'LineWidth',2, 'Color', 'k');

decay_start = 22;
plot(decay_start:n, exp_decay(decay_start,2:(n-decay_start + 2)), '-', 'LineWidth',2, 'Color', 'k');

% Loss model
% Specify the Loss Rate
% newExp = max(min_exp, currentExp-loss_rate*(1 - fbValue)*alpha*(1-currentExp));
loss_rate = 2;

exp_loss = zeros(n, n);
exp_loss(:, 1) = exp_increase;
for i=2:n
    for j=1:n-1
        exp_loss(i, j+1) = max(min_exp, exp_loss(i, j)-loss_rate*(1-uncoop_interaction(j))*alpha*(1-exp_loss(i, j)));
    end
end
loss_start = 5;
exp_temp = exp_loss(loss_start, (exp_loss(loss_start, :) > 0));
exp_temp = [exp_temp 0];
plot(loss_start:(loss_start+size(exp_temp, 2)-1), exp_temp, '-hr', 'LineWidth',1, 'MarkerEdgeColor','b', 'MarkerFaceColor','w');

loss_start = 22;
exp_temp = exp_loss(loss_start, (exp_loss(loss_start, :) > 0));
exp_temp = [exp_temp 0];
plot(loss_start:(loss_start+size(exp_temp, 2)-1), exp_temp, '-sr', 'LineWidth',1, 'MarkerEdgeColor','b', 'MarkerFaceColor','w');

xlabel('Time', 'FontSize',13);
ylabel('Experience Value', 'FontSize',13);
title('Experience Model consists of Develop, Decay and Loss trends', 'FontSize',15);
legend({'Develop','Decay (Weak Tie)', 'Decay (Strong Tie)', 'Loss (Weak Tie)', 'Loss (Strong Tie)'},'FontSize',15);