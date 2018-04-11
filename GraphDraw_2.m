% Graph Draw for the paperwork
figure('NumberTitle', 'on', 'Name', 'This is the figure title');
x = linspace(1, 100, 100);
title('Comparison between Trust-based, Previous-Best-QoD-based and Average-QoD-based MCS Schemes', 'FontSize',15);


subplot(2,2,1)       % add first plot in 2 x 2 grid
plot(x, -result_100tasks/100, '--*', 'MarkerSize',8);
grid on;
grid minor;
xlabel('Number of Malicious Users', 'FontSize',12);
ylabel('Damage Level', 'FontSize',12);         
title('Number of Tasks: 100')
legend({'Trust-based Scheme','PrevBest-QoD-based Scheme', 'Average-QoD-based Scheme'},'FontSize', 10, 'Location','northwest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,2)       % add second plot in 2 x 2 grid
plot(x, -result_400tasks/400, '--*', 'MarkerSize',8);  
grid on;
grid minor;
xlabel('Number of Malicious Users', 'FontSize',12);
ylabel('Damage Level', 'FontSize',12);         
title('Number of Tasks: 400')
legend({'Trust-based Scheme','PrevBest-QoD-based Scheme', 'Average-QoD-based Scheme'},'FontSize', 10, 'Location','northwest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3)       % add third plot in 2 x 2 grid
plot(x, -result_800tasks/800, '--*', 'MarkerSize',8);           
grid on;
grid minor;
xlabel('Number of Malicious Users', 'FontSize',12);
ylabel('Damage Level', 'FontSize',12);         
title('Number of Tasks: 800')
legend({'Trust-based Scheme','PrevBest-QoD-based Scheme', 'Average-QoD-based Scheme'},'FontSize', 10, 'Location','northwest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,4)       % add fourth plot in 2 x 2 grid
plot(x, -result_1600tasks/1600, '--*', 'MarkerSize',8);  
grid on;
grid minor;
xlabel('Number of Malicious Users', 'FontSize',12);
ylabel('Damage Level', 'FontSize',12);         
title('Number of Tasks: 1600')
legend({'Trust-based Scheme','PrevBest-QoD-based Scheme', 'Average-QoD-based Scheme'},'FontSize', 10, 'Location','northwest');

h = suptitle('Comparisons among Trust-based, Previous-Best-QoD-based, and Average-QoD-based Schemes');
set(h,'FontSize',14,'FontWeight','bold');
