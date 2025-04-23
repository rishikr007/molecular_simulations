function [pi_approx] = monte_carlo_pi(n_trials)
    x = rand(n_trials, 1);
    y = rand(n_trials, 1);
    n_hits = sum((x.^2 + y.^2) <= 1);
    pi_approx = 4 * n_hits / n_trials;
endfunction

num_trials = 500;
pi_values = zeros(num_trials, 1);

for i = 1:num_trials
    pi_values(i) = monte_carlo_pi(i);
end

//disp(pi_values)

plot(1:num_trials, pi_values, 'b-', 'LineWidth', 2);
xlabel("Number of Trials");
ylabel("Approximation of π");
title("Monte Carlo Approximation of π");
