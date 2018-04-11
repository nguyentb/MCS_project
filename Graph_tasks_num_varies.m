x = linspace(1, 160, 160);
y = linspace(0, 160, 9);
grid on;
hold on;
grid minor;
xlabel('Number of Requested Services', 'FontSize',12);
ylabel('QoS score', 'FontSize',12);

plot(x, result_draw(:,1), '--*', 'MarkerSize',8);
plot(x, result_draw(:,3), '--*', 'MarkerSize',8);
plot(x, result_draw(:,4), '--*', 'MarkerSize',8);
plot(x, result_draw(:,2), '--*', 'MarkerSize',8);
%plot(x, result_draw2(:,5), '--*', 'MarkerSize',8);

xlim([0 160]);
ylim([1.5 4]);
xticks(y);
xticklabels({y});

%title('QoS vs Number of Services in different User Recruitment schemes', 'FontSize',14);
legend({'Trust-based Scheme', 'Average-QoD-based Scheme', '3-degree Polynomial Regression', 'Random Selection Scheme'},'FontSize', 9, 'Location','northwest', 'FontSize',10);