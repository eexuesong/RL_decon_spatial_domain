%% Function: 1-D Guaussian fitting
function [sigma, mu, A] = mygaussfit(x, y, h)
% [sigma,mu,A] = mygaussfit(x,y)
% [sigma,mu,A] = mygaussfit(x,y,h)

%{
This function is doing fit to the function
y = A * exp(-(x - mu)^2 / (2 * sigma^2))

The fitting is been done by a polyfit the lan of the data.

h is the threshold which is the fraction from the maximum y height that the data is been taken from.
h should be a number between 0-1.
If h have not been taken it is set to be 0.2 as default.
%}

% threshold
if nargin == 2
    h = 0.2;
end

% cutting
ymax = max(y);
xnew = [];
ynew = [];
for n = 1:length(x)
    if y(n) > ymax * h
        xnew = [xnew, x(n)];
        ynew = [ynew, y(n)];
    end
end

% fitting
ylog = log(ynew);
xlog = xnew;
p = polyfit(xlog, ylog, 2);
A2 = p(1);
A1 = p(2);
A0 = p(3);

%{
Note: ln(y) = -1 / (2 * sigma ^2) * x^2 + mu / sigma^2 * x + ln(A) - mu^2 / (2 * sigma^2)
% A2 = -1 / (2 * sigma ^2)
% A1 = mu / sigma^2
% A0 = ln(A) - mu^2 / (2 * sigma^2)
%}

sigma = sqrt(-1/(2 * A2));
mu = A1 * sigma^2;
A = exp(A0 + mu^2 / (2 * sigma^2));
end