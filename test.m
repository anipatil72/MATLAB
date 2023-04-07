% Define the first function to plot
f1 = @(x) sin(x);

% Create a new figure window for the first plot
figure(1)

% Plot the first function in a single subplot
subplot(1, 1, 1)
fplot(f1, [-pi, pi])
title('Plot of sin(x)')

% Define the second function to plot
f2 = @(x) x.^2;

% Create a new figure window for the second plot
figure(2)

% Plot the second function in a single subplot
subplot(1, 1, 1)
fplot(f2, [-5, 5])
title('Plot of x^2')
