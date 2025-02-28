function dydx = FiniteDiff(y,h,varagin)

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
