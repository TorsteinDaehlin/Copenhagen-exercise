function dydx = FiniteDiff(y,h,varagin)
% finiteDiff.m
% -------------------------------------------------------------------------
% Differentiates the equally spaced dependent variables y = f(x) using a
% finite difference scheme. The function can return the first or second
% derivative of y.
% -------------------------------------------------------------------------
% Syntax and description:
% dydx = finiteDiff(y,h) returns the first derivative of the input array y
% = f(x) with equally spaced steps specified by the step size h.
%
% dydx = finiteDiff(y,h,order) returns the derivative of the specified
% order. The input argument 'order' can take values 1 or 2, returning the
% first or second order derivative, respectively.
%
% The first-order method utilizes a two-sided two-point scheme to calculate
% intermediate points, while one-sided forward and backward two-point
% schemes is used to calculate the first and last data point, respectively.
%
% The second-order method utilizes a two-sided three-point scheme to
% calculate intermediate points, while one-sided forward and backward
% three-point schemes is used to calculate the first and last data point,
% respectively.
% -------------------------------------------------------------------------
% Written by Torstein E. Daehlin, August 2019
% -------------------------------------------------------------------------

% Detect and assigne variable input arguments
nargs = nargin;
if nargs == 3
    order = varagin(1);
else
    order = 1;
end

% Determine size of input
[row, col] = size(y);

% Transpose input if there are more columns than rows
transpose = false;
if col > row
    y = y';
    [row, col] = size(y);
    transpose = true;
end

% Preallocate output array
dydx = zeros(row, col);

% Select first or second order case
if order == 1
    for cdx = 1:col
        % Differentiate first point using one-sided forward two-point scheme
        dydx(1,cdx) = (y(2,cdx) - y(1,cdx))/h;

        % Differentiate 2:n-1 points using two-sided two-point scheme
        for idx = 2:(row-1)
            dydx(idx,cdx) = (y(idx+1,cdx) - y(idx-1,cdx))/(2*h);
        end

        % Differentiate point n using one-sided backward two-point scheme
        dydx(idx+1,cdx) = (y(idx+1,cdx) - y(idx,cdx))/h;
    end
elseif order == 2
    for cdx = 1:col
        % Differentiate first point using one-sided forward three-point scheme
        dydx(1,cdx) = (y(3,cdx) - 2*y(2,cdx) + y(1,cdx))/h^2;

        % Differentiate 2:n-1 points using two-sided three-point scheme
        for idx = 2:(row-1)
            dydx(idx,cdx) = (y(idx+1,cdx) - 2*y(idx,cdx) + y(idx-1,cdx))/h^2;
        end

        % Differentiate point n using one-sided backward two-point scheme
        dydx(idx+1,cdx) = (y(idx+1,cdx) - 2*y(idx,cdx) + y(idx-1,cdx))/h^2;
    end
else
    error('The input argument ''order'' must take values 1 or 2');
end

% If input data series had time in columns, transpose the output
if transpose
    dydx = dydx';
end
end
