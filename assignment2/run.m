%% problem1
R = 20; %repeats for system 1
ap = zeros(1,R); %holds average price for system 1
ad = zeros(1,R); %holds average delay for system 1

for i = 1:R
    [ap(i), ad(i), ~] = shopSimulation(5,true);
end

%according to Bonferroni inequality, in order to be at least 1 - a confidence interval
 %each separate interval should be 1 - a/c, where c is the number of
 %measures to make
 a_overal = 0.5;
 c = 2; %the number of estimators 
 a_ind = a_overal/c;
ci_ap = get_ci(ap, a_ind);
ci_ad = get_ci(ad, a_ind);
fprintf("The 95%% overall confidence interval for average profit: [%.2f, %.2f], and for average delay: [%.2f, %.2f]\n", ...
    ci_ap(1), ci_ap(2), ci_ad(1), ci_ad(2))
% null hypothesis: the mean of average profit is 300
% For a single run, we got the [291.50, 300.50] for the average profit interval, this interval includes 300, then
%  with significance level of 5%, we cannot conclude the mean is different from 300
% 


%% problem2
R2 = 30; %the repeats of system 2
ap2 = zeros(1,R2); %holds the average price for system2
a2 = 0.05; %alpha  = 1 - confidence interval
for i = 1:R2
    [ap2(i), ~, ~] = shopSimulation(5,false); %set the second parameter to false for system 2
end

zbar = mean(ap) - mean(ap2); %the difference in mean as point estimate
vara = var(ap)/length(ap); %the variance for the mean of system 1
varb = var(ap2)/length(ap2); %the variance for the mean of system 2
f = round((vara + varb)^2/(vara^2/(R-1) + varb^2/(R2-1))); % the degree of freedom for the unpaired t test
ci = [zbar + tinv(a2/2,f)*sqrt(vara+varb), zbar + tinv(1-a2/2,f)*sqrt(vara+varb)]; %one single run of 95% confidence interval 
fprintf("95%% confidence interval: [%.2f %.2f]\n", ci(1), ci(2));

t = zbar/sqrt(vara + varb); %the t-statistic
p = 2 *(1 - tcdf(t, f)); %the p value, needs to double for both ends
fprintf("p value: %.2f%%\n",p*100);

%h0 the null hypothesis is there is no difference between the two methods,
%i.e. u1 = u2
% the t-statistic for the single run is 3.2955
% the 95% confidence interval run for the single run is [5.25 22.51], this
% does not include the 0 value, so there is difference in the means of of 2 systems. Similarly, the p value for the single run
% is 0.27% which is less than 5%, and the same conclusion can be drawn.

%% problem 3
%running takes time, please be patient
capacity = [5 5 7 3];
strategy = [true false true true];
k = length(capacity); %the number of configurations
N = zeros(k,1); %store the total number of obversation for each config
Ssquared = zeros(k,1); %store the sample variances of the initial run
Xbar1 = zeros(k,1); %store average of the initial run for each config
Xtilde = zeros(k,1); %store the weighted sample mean for each config
n0 = 20; %the initial fixed run
d = 0.4; %the distance 
h1 = 1.896; %found at the appendix with p = 0.9 and n0 = 20
for i = 1:k
    x = zeros(n0, 1); %store the current averages
    for j = 1:n0
        [x(j), ~, ~] = shopSimulation(capacity(i), strategy(i));
    end
    varx = var(x); %var of average
    Ssquared(i) = varx;
    Xbar1(i) = mean(x);
    n_total = max(n0+1,ceil(h1^2*varx/(d^2))); %total runs needed for current config
    N(i) = n_total;
    w1 = n0/n_total*(1 + sqrt(1 - n_total/n0 * (1 - (n_total - n0) * d^2/(h1^2*var(x))))); %weight 1 for phase 1 mean
    w2 = 1 - w1; %weight 2 for phase 2 mean
    x2 = zeros(n_total - n0, 1); %store the averages for the phase 2 mean
    for k = 1:n_total - n0
        [x2(k), ~, ~] = shopSimulation(capacity(i), strategy(i));
    end
    Xtilde(i) = w1 * mean(x) + w2 * mean(x2);    
end

Xtilde

% result of one run:
%   286.3475
%   279.9684
%   225.5851
%   313.3167
% based on that, we choose strategy 4 since this is the most profitable.


%%
function ci = get_ci(x, a)
    se = std(x)/sqrt(length(x)); %standard error               
    ts = tinv([a/2  1-a/2],length(x)-1);     %lo=a/2, hi=1-a/2, df=n-1
    ci = mean(x) + ts*se; 
end
    

