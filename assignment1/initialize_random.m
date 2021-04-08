function[f] = initialize_random(seed, a, c, m)
%initialize lcg generator
    x0 = seed;
    function [val] = g()
         val = mod(a*x0 + c, m);
         x0 = val;
    end
    f = @g;
end