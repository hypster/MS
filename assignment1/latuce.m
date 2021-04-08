generate = initialize_random(100,18,0,101);
generate2 = initialize_random(100,2,0,101);

% simulate(generate);
% simulate(generate2);

seed = 30;
a = 25214903917;
c = 11;
m=2^24;
l = 10000;
alpha = 1.358 % alpha for 95%confidence interval of known distribution

% 
% out = zeros(1,l);
% generate3 = initialize_random(seed, a, c, m);
% simulate(generate3,[]);

% generate4 = initialize_random(100,65539,0,2^31);
% simulate(generate4,[]);

o = simulate(initialize_random(1,7,2,11), []);


function [o] = simulate(generate, n)
o = zeros(1,2000);
for i = 1:2000
    o(i) = generate();
end

if n == 3
    o3 = o(3:end);
    o2 = o(2:end-1);
    o1 = o(1:end-2);
    figure;
    scatter3(o1,o2,o3);
else 
    o2 = o(2:end);
    o1 = o(1:end-1);
    figure;
    scatter(o1,o2);
end
end


