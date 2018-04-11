% Graph Draw for the paperwork
figure('NumberTitle', 'on', 'Name', 'This is the figure title');
title('Comparison between different User Recruitment Schemes in MCS', 'FontSize',15);
x = linspace(1, 100, 100);
y = linspace(0, 25, 6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1)       % add first plot in 2 x 2 grid
plot(x, result_draw1, '--*', 'MarkerSize',8);
xticks([0 20 40 60 80 100]);
xticklabels({'0%', '5%', '10%', '15%', '20%', '25%'});
ylim([1.5 4]);

grid on;
grid minor;
xlabel('Percentage of Malicious Users', 'FontSize',12);
ylabel('QoS Score', 'FontSize',12);         
title('Number of Requested Services: 10')
legend({'Trust-based', 'Average-QoD-based', '3-degree Polynomial Regression', 'Random Selection'},'FontSize', 9, 'Location','southeast');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,2)       % add second plot in 2 x 2 grid
plot(x, result_draw2, '--*', 'MarkerSize',8);  
xticks([0 20 40 60 80 100]);
xticklabels({'0%', '5%', '10%', '15%', '20%', '25%'});
ylim([1.5 4]);

grid on;
grid minor;
xlabel('Percentage of Malicious Users', 'FontSize',12);
ylabel('QoS Score', 'FontSize',12);         
title('Number of Requested Services: 40')
legend({'Trust-based', 'Average-QoD-based', '3-degree Polynomial Regression', 'Random Selection'},'FontSize', 9, 'Location','southeast');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3)       % add third plot in 2 x 2 grid
plot(x, result_draw3, '--*', 'MarkerSize',8);
xticks([0 20 40 60 80 100]);
xticklabels({'0%', '5%', '10%', '15%', '20%', '25%'});
ylim([1.5 4]);

grid on;
grid minor;
xlabel('Percentage of Malicious Users', 'FontSize',12);
ylabel('QoS Score', 'FontSize',12);         
title('Number of Requested Services: 80')
legend({'Trust-based', 'Average-QoD-based', '3-degree Polynomial Regression', 'Random Selection'},'FontSize', 9, 'Location','southeast');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,4)       % add fourth plot in 2 x 2 grid
plot(x, result_draw4, '--*', 'MarkerSize',8);
xticks([0 20 40 60 80 100]);
xticklabels({'0%', '5%', '10%', '15%', '20%', '25%'});
ylim([1.5 4]);

grid on;
grid minor;
xlabel('Percentage of Malicious Users', 'FontSize',12);
ylabel('QoS Score', 'FontSize',12);         
title('Number of Requested Services: 160')
legend({'Trust-based', 'Average-QoD-based', '3-degree Polynomial Regression', 'Random Selection'},'FontSize', 9, 'Location','southeast');

%h = suptitle('Changes in Percentage of Malicious Users in different User Recruitment Schemes');
%set(h,'FontSize',16,'FontWeight','bold');