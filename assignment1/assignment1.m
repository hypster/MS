%% starter
close all;clear all;
load('dataAss1.mat');
%% data analysis
figure;
boxplot(dat0);
% the one outliner is too extreme, remove it
[~,i] = max(dat0);
dat0(i) = [];
figure;
qqplot(dat0); 
figure;
hist(dat0,20);

%summary stats
n = length(dat0);
u = mean(dat0);
sigma = std(dat0);
m = median(dat0);
sk = n^2/(n-1)/(n-2) * sum((dat0-u).^3)/n/sigma^3; %the skewness is corrected for bias
fprintf("mean = %.4f, standard deviation = %.4f, \n median = %.4f, skewness = %.4f\n", u,sigma,m, sk);

%from the plot, we can see that the distibution is not symmetric, here I choose weibull distribution, since it fits good to a skewed distribution with tail to the right and is versatile
%% Parameter estimaton using MLE reference textbook 290-292
n = length(dat0);
ak = ((6/pi^2) * (sum(log(dat0).^2) - (sum(log(dat0)))^2/n) / (n-1))^(-1/2);
A = sum(log(dat0))/n;

epi = 1e-6;
cnt = 0;
while (1)
    cnt = cnt+1;
    old = ak;
    Bk = sum(dat0.^ak);
    Ck = sum(dat0.^ak .* log(dat0));
    Hk = sum(dat0.^ak .* power(log(dat0),2));
    ak = ak + (A + 1/ak - Ck/Bk)/(1/ak^2 + (Bk*Hk-Ck^2)/Bk^2);
    if abs(ak-old) < epi
        fprintf("iteration = %d\n", cnt);
        break
    end
end
beta = (sum(dat0.^ak)/n)^(1/ak);
alpha = ak;
fprintf("alpha = %.4f, beta = %.4f\n", ak, beta)
%% density histogram
fx = @(x) alpha *beta^(-alpha) * x.^(alpha-1).*exp(-(x/beta).^alpha);
xs = min(dat0):0.1:max(dat0);
figure;
histogram(dat0, 'Normalization','pdf');
hold on;
% n = length(dat0);
plot(xs, fx(xs),'linewidth',5);
title("density histogram plot");
% as we can see from the graph, the hypothesis seems ok
%% chi-squared test
%h0: the data has a weibull distribution
%the bins are not combined for size less than 5, if use toolbox function set Emin = 0 will get
%similar value
n = 53; %the number of bins chosen, this number satisfies the recommendation range [sqrt(n), n/5] for n > 100
p = 1/n;
gi = @(i) (-log(1- i*p))^(1/alpha) * beta; % this calculates the endpoint for equal prob interval
itvl  = zeros(1, n);
itvl(n) = inf;
for i = 1:n-1
    itvl(i) = gi(i);
end
dat0_s = sort(dat0);
fq = zeros(1,n);
pos = 1;
for i = 1: length(dat0)
    if dat0_s(i) < itvl(pos)
        fq(pos) = fq(pos)+1;
    else
        while itvl(pos) <= dat0_s(i)
            pos = pos+1;
        end
        fq(pos) = fq(pos)+1;
    end
end

ex = length(dat0) /n;
chisq = sum((fq - ex).^2)/ex;
a_50_05 = 67.5;
if chisq < a_50_05
    fprintf("chi-square  = %.2f, less than alpha(50,0.05) = %.2f, h0 is not rejected\n", chisq, a_50_05);
else 
    fprintf("chi-square  = %.2f, greater than alpha(50,0.05) = %.2f, h0 is rejected\n", chisq, a_50_05);
end
%% optional, see each interval represents equal probability
% Fx  = @(x) 1- exp(-(x/beta).^alpha);
% test = zeros(1, length(itvl));
% for i = 1:length(itvl)
%     edpt = itvl(i)
%     test(i) = Fx(edpt);
% end
% for i = length(test):-1:2
%     test(i) = test(i) -test(i-1);
% end
% plot(test)
%% truncated lcg
seed = 999999;
a = 25214903917;
c = 11;
m=2^48;
l = 10000;
lo = 15;
hi = 47;
alpha = 1.358; % alpha for 95% confidence interval of known distribution

out = zeros(1,l);
generate = initialize_random(seed, a, c, m, lo, hi);
for i = 1:l
    out(i) = generate();
%     out(i) = rand();
end


%% poker
E = [0.3024 0.5040 0.1080 0.072 0.009 0.0045 0.0001]; %probability for the 7 patterns

cnt = zeros(1,7); % store the frequency for each pattern

for i = 1:5:length(out)
    lookup = zeros(10,5); %table holds positions for all possible samples in a hand of 5 cards experiment
    hand = out(i:i+4); %get the next 5 cards
    row = floor(hand * 10) + 1; %transform the random number into corresponding positions in the table
    col = 1:5;
    ind = 10*(col-1) + row; %use here matlab unique index method
    lookup(ind) = 1; %set corresponding positions to occupied
    bin_num = sum(lookup,2); %get row sum, which represents the number for each type of card in one hand
    if sum(bin_num == 5)==1 %five of a kind
        cnt(7) = cnt(7)+1;
    elseif sum(bin_num == 4) == 1 %four of a kind
        cnt(6) = cnt(6)+1;
    elseif sum(bin_num == 3) == 1 
        if sum(bin_num==2) == 1 %full house
            cnt(5) = cnt(5)+1;
        else 
            cnt(4) = cnt(4)+1; %three of a kind
        end
    elseif sum(bin_num == 2) == 2
        cnt(3) = cnt(3)+1; %two pairs
    elseif sum(bin_num==2) == 1
        cnt(2) = cnt(2)+1; %one pair
    else
        cnt(1) = cnt(1)+1; %all different
    end
end

E = E * 2000; %np
X2 = sum((cnt-E).^2 ./E);
alpha_df6_05 = 12.59; %df = number of categories - 1
if(X2 < alpha_df6_05)
    fprintf('chi2 = %.2f, less than a(6,0.05) = %.2f, null hypothesis can not be rejected.\n', X2, alpha_df6_05);
else
    fprintf('chi2 = %.2f, greater than a(6.0.05) = %.2f, null hypothesis is rejected. The generated number fails the independence test\n', X2, alpha_df6_05);
end

    
%% ks test
out = sort(out);

% h0: the generated number is sampled from a uniform distribution
x = (1:l)/l; %i/n
x2 = (0:l-1)/l; %(i-1)/n
d_plus = max(x - out);
d_min = max(out - x2);
d = max(d_plus,d_min);
adjusted_d = (sqrt(l) + 0.12 + 0.11/sqrt(l)) * d; 
figure;
plot([0,1], [0,1], 'r--');
hold on;
stairs(out, x, 'b'); %empirical distribution
title("comparison between uniform and empirical distribution");
%we can see that there is little difference between the two lines
if(adjusted_d < alpha)
    fprintf('d = %.2f, less than a = %.2f, null hypothesis can not be rejected.\n', adjusted_d, alpha);
else
    fprintf('d = %.2f, greater than a = %.2f, null hypothesis is rejected. The generator does not satisfy uniform distribution\n', adjusted_d, alpha);
end


%%
function[f] = initialize_random(seed, a, c, m, lo, hi)
%return a function for generating random number with parameters
% the generated integer is between 0 and 2^(hi-lo) paper reference: https://doi.org/10.1007/BFb0052242
    x0 = seed;
    function [k] = g()
         val = mod(a*x0 + c, m);
         x0 = val;
         d = mod(val, 2^lo); %truncate lower bits
         k = mod(val - d, 2^hi)/(2^lo); %k is an integer because of the previous mod operation
         k = k/2^(hi-lo); % return the random number in u(0,1)
    end
    f = @g;
end

