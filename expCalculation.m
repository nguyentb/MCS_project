function [ newExp ] = expCalculation( fbValue, currentExp, previousExp)
% Experience Calculation

%   Supportive Feedback: fbValue >= 0.85
%   Neutral Feedback: 0.75 <= fbValue < 0.85
%   Unsupportive Feedback: fbValue < 0.75

    supportive_threshold = 0.7;
    unsupportive_threshold = 0.5;

    % Development Model: let define alpha = maximum increase of Experience
    % x(i+1)= (1 - alpha)*x(i) + alpha;
    alpha = 0.15;
    if (fbValue >= supportive_threshold)
        newExp = currentExp + fbValue*alpha*(1-currentExp);
    end
    
    % Decay Model: let define alpha = maximum increase of Experience
    % exp(t+1) = exp(t) - decay
    min_decay = 0.005;
    decay_rate = 0.005;
    exp_bootstrap = 0.5;        
    if (fbValue < supportive_threshold) && ((fbValue >= unsupportive_threshold))        
        newExp = max(exp_bootstrap, currentExp - min_decay*(1 + decay_rate - previousExp));        
    end

    % Loss model
    % Specify the Loss Rate
    min_exp = 0;
    loss_rate = 3;
    if (fbValue < unsupportive_threshold)        
        newExp = max(min_exp, currentExp-loss_rate*(1 - fbValue)*alpha*(1-currentExp));
    end
end