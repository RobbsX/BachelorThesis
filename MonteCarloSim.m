function [slope, probPosRes] = MonteCarloSim(tradesAllClean, plott)


data = tradesAllClean;

[f,x] = ecdf(data);
Tempx = abs(x-0);
closesttozero = Tempx == min(Tempx);
probPosRes = 1-f(closesttozero);
probPosRes = probPosRes(end); % take the worst possibility

ECnew = emprand(data, length(data), 1000); 
ECnew = [zeros(1, 1000); ECnew];
ECnew = cumsum(ECnew);
qu5 = quantile(ECnew', (0.05));
qu95 = quantile(ECnew', (0.95));
avg = quantile(ECnew', (0.5));

if plott == 1
    figure; histogram(data, 50);
    figure; ecdf(data);
    figure;plot(ECnew);
    hold on;
    plot(data, 'b-','LineWidth', 3);
    hold on;
    plot(qu5, 'LineWidth', 5);
    hold on;
    plot(avg, 'LineWidth', 5);
    hold on;
    plot(qu95, 'LineWidth', 5);
end

% Find out Slope of average using Linear Regression
linLen = (1:length(avg));
[p] = polyfit(linLen,avg, 1);
slope = p(1);

% Je risikoaverser desto höher die untere Quartile (z.B. 10%)
