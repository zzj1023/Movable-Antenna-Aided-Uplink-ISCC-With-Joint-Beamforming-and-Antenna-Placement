clc
clear
% 鍔犺浇鏁版嵁锛堜繚鎸佷笉鍙橈級
totalsum_proposed = load(strcat('..\res_P\proposed_SDR_comput.mat'), 'totalsum_proposed');
totalsum_proposed = totalsum_proposed.totalsum_proposed;
totalsum_GA = load(strcat('..\res_P\GA_SDR_2.mat'), 'GA_sumRate');
totalsum_GA = totalsum_GA.GA_sumRate;
totalsum_FPA = load(strcat('..\res_P\FPA_SDR_2.mat'), 'FPA_sumRate');
totalsum_FPA = totalsum_FPA.FPA_sumRate;
% totalDelay_onlyMA = load(strcat('..\Iter_50\onlyMA.mat'), 'totalDelay_onlyMA');
% totalDelay_onlyMA = totalDelay_onlyMA.totalDelay_onlyMA;
totalsum_AS = load(strcat('..\res_P\AS_SDR_2.mat'), 'finalRates_AS');
totalsum_AS = totalsum_AS.finalRates_AS;
totalsum_Random = load(strcat('..\res_P\Random_SDR_2.mat'), 'Random_sumRate');
totalsum_Random = totalsum_Random.Random_sumRate;

figure;
x = [0.2 0.4 0.6 0.8 1 1.2].';

set(gcf, 'DefaultTextFontName', 'Times New Roman', 'DefaultTextInterpreter', 'tex'); 

plot(x,totalsum_proposed, '-k', 'LineWidth', 1);
h_proposed = plot(x,totalsum_proposed, '-ok', 'MarkerFaceColor', 'w', 'LineWidth',1); % 闂撮殧鏍囪
hold on;

plot(x,totalsum_GA, '--g', 'LineWidth', 1);
h_GA = plot(x,totalsum_GA, '>g', 'MarkerFaceColor', 'w', 'LineWidth',1);

plot(x,totalsum_FPA, '--r', 'LineWidth', 1);
h_FPA = plot(x,totalsum_FPA, '^r', 'MarkerFaceColor', 'w', 'LineWidth',1);


plot(x,totalsum_AS, '--b', 'LineWidth', 1);
h_AS = plot(x,totalsum_AS, 'db', 'MarkerFaceColor', 'w', 'LineWidth',1);

plot(x,totalsum_Random, '--m', 'LineWidth', 1);
h_Random = plot(x,totalsum_Random, 'sm', 'MarkerFaceColor', 'w', 'LineWidth',1);

% 鍥句緥锛堜繚鎸佷笉鍙橈級
% h_legend = legend([h_proposed,h_FPA,h_AS], 'Proposed','FPA','AS', 'Location','southeast');
h_legend = legend([h_proposed,h_GA,h_FPA,h_AS,h_Random], ...
    'Proposed','GA','FPA','AS','RPA', 'Location','southeast');
% h_legend = legend([h_proposed,h_FPA,h_onlyMA,h_AS], 'Proposed','FPA','PSO-MRT','AS', 'Location', 'southeast');
set(h_legend, 'FontSize', 12);

xlabel('The maximum transmit power (W)', 'FontSize', 14);
ylabel('Sum rate (bps/Hz)', 'FontSize', 14);
xticks([0.2 0.4 0.6 0.8 1 1.2]);
xlim([0.2,1.2]);
grid on;
ax1 = gca;
ax1.GridLineStyle = '--';
