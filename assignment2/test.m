a = [2.1270 4.2164 4.3521 5.4543 3.5270];
b = [4.5409 4.7218 4.1447 1.8607 1.3491];
z = a-b;

% [-1.6196, 2.8434]
zbar = mean(z)
% varz = sum((z - zbar).^2)/(5*4)
varz = var(z)/5

[zbar + tinv(0.05,4)*sqrt(varz), zbar + tinv(0.95,4)*sqrt(varz)]

%%
b = [b 2.0290 8.9853 7.4060 7.2067 4.2002]
%%
% welch [-2.4396, 1.0214]
abar = mean(a);
bbar = mean(b);
zbar = abar - bbar;
vara = var(a)/length(a);
varb = var(b)/length(b);
f = (vara + varb)^2/(vara^2/(length(a)-1) + varb^2/(length(b)-1))
[zbar + tinv(0.05,f)*sqrt(vara+varb), zbar + tinv(0.95,f)*sqrt(vara+varb)]

%%
% \tilde{X}z = 4.0365 \tilde{X}k=3.7553 Since we want to minimise the delay we choose Klunkytel
close all; clear; 
X = [2.1270 4.5409;
4.2164 4.7218;
4.3521 4.1447;
5.4543 1.8607;
3.5270 1.3491;
1.4878 2.0290;
1.5965 8.9853;
5.3555 7.4060;
5.0970 7.2067;
1.8890 4.2002;
1.2391 2.6969;
6.2157 1.1058;
5.2080 10.3219;
1.8122 0.7468;
10.7755 4.3590;
1.2525 11.9941;
6.2235 1.6067;
2.0080 8.3922;
2.7263 3.0037;
3.8331 2.3109];


n0 = size(X,1);
ns = zeros(1,2);
h1 = 1.896;
d = 0.4;
x2 = [4.0844  3.6597];
res = zeros(1,2);
for i = 1:size(X,2)
    x = X(:,i);
    n0 = length(x);
    varx = var(x);
    n_total = max(n0+1,ceil(h1^2*varx/(d^2)));
    w1 = n0/n_total*(1 + sqrt(1 - n_total/n0 * (1 - (n_total - n0) * d^2/(h1^2*var(x)))));
    w2 = 1 - w1;
    res(i) = w1*mean(x) + w2 * x2(i);
end

res

