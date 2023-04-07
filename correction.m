% Define the size of the pupil grid
n = 256;
m = n;

% Create a grid of x and y coordinates
x = linspace(-1, 1, n);
y = linspace(-1, 1, m);
[X, Y] = meshgrid(x, y);
R = sqrt(X.^2 + Y.^2);
Theta = atan2(Y, X);

% Generate a pupil function with tilt
tilt = 0.05;
pupil = zeros(n, m);
for i = 1:n
    for j = 1:m
        pupil(i, j) = exp(-((X(i, j)-tilt*Y(i, j))^2 + Y(i, j)^2)/(2*(0.1)^2));
    end
end

% Define the Zernike order and frequency to correct for tilt
n = 2;
m = 2;

% Compute the Zernike polynomial for tilt correction
Z = zernfun2(n, m, R, Theta);

% Compute the coefficients of the Zernike polynomial
c = Z(:)\pupil(:);

% Reconstruct the phase aberration using the Zernike coefficients
pupil_corrected = reshape(Z*c, [size(pupil, 1), size(pupil, 2)]);

% Plot the original and corrected pupil functions
figure;
subplot(1,2,1);
imagesc(pupil); axis square; colorbar;
title('Original Pupil Function');
subplot(1,2,2);
imagesc(real(pupil_corrected)); axis square; colorbar;
title('Corrected Pupil Function');



function Z = zernfun2(n, m, R, Theta)
% Compute the Zernike polynomial with given order and frequency at each point in the grid
% Inputs:
%   n: order of the Zernike polynomial
%   m: frequency of the Zernike polynomial
%   R: matrix of radial distances for each point on the grid
%   Theta: matrix of azimuthal angles for each point on the grid
% Output:
%   Z: matrix of Zernike polynomial values for each point on the grid

% Check that the inputs are valid
if (n < 0) || (mod(n-m, 2) ~= 0)
    error('Invalid parameters for Zernike polynomial.');
end

% Compute the normalization constant for the Zernike polynomial
norm_factor = sqrt((2*(n+1))/(1+double(m==0)));

% Initialize the Zernike polynomial to zero
Z = zeros(size(R));

% Compute the Zernike polynomial for each value of R and Theta
for k = 0:((n-m)/2)
    num = (-1)^k * factorial(n-k);
    denom = factorial(k) * factorial((n+m)/2-k) * factorial((n-m)/2-k);
    term = num / denom * R.^(n-2*k);
    Z = Z + term .* exp(1i*m*Theta);
end

% Apply the normalization factor to the Zernike polynomial
Z = norm_factor * Z;
end

