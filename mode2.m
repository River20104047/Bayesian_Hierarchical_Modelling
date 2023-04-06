function [mod] = mode2(dat)
%   Calculation of mode based on kernel density
% Kernel density
TFNaN = sum(isnan(dat(:))); % Number of NaN
TFInf = sum(isinf(dat(:)));% Number of Inf
if TFNaN + TFInf > 0
    mod = NaN;
else   
    a      = prctile(dat,15);
    b      = prctile(dat,85);
    c      = (b - a) / 1000;
    pts    = a:c:b;
    [f,xi] = ksdensity(dat,pts);
    % plot(xi,f) % plot kernel density if needed
    % Find xi that gives max f
    m      = max(f);
    in     = find(f == m);
    mod    = xi(in);
end

