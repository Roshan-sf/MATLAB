% MATLAB Script to Create Subplots for Railgun Power Source Trade Study

% Data from trade study table
power_sources = {'Battery Bank', 'Capacitor Bank', 'Hybrid System'};
energy_efficiency = [85, 95, 90];  % Percentage
power_density = [1.2, 5.5, 3.5];    % kW/kg
recharge_time = [45, 5, 15];        % Minutes
cost_per_kw = [500, 1200, 800];     % $ per kW

% Create figure with subplots
figure;

% Subplot 1: Energy Efficiency
subplot(2,2,1);
bar(categorical(power_sources), energy_efficiency);
ylabel('Efficiency (%)');
title('Energy Efficiency Comparison');
grid on;

% Subplot 2: Power Density
subplot(2,2,2);
bar(categorical(power_sources), power_density);
ylabel('Power Density (kW/kg)');
title('Power Density Comparison');
grid on;

% Subplot 3: Recharge Time
subplot(2,2,3);
bar(categorical(power_sources), recharge_time);
ylabel('Recharge Time (min)');
title('Recharge Time Comparison');
grid on;

% Subplot 4: Cost per kW
subplot(2,2,4);
bar(categorical(power_sources), cost_per_kw);
ylabel('Cost ($/kW)');
title('Cost Comparison');
grid on;

% Adjust layout for better visibility
sgtitle('Railgun Power Source Trade Study');
