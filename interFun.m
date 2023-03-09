function [Zs, dZs] = interFun(x, limits)
%% Interpolation scheme based on power-law approach
% 
% Inputs:
%   x       - design variable [0,1]^Tx1
%   limits  - boundary values for the interpolated resistivities
% 
% Outputs:
%   Zs      - interpolated surface resistivities
%   dZs     - derivative of the interpolation function  
% 
% 2023, Jonas Tucek, CTU in Prague, jonas.tucek@fel.cvut.cz

Zvac = limits(1);
Zmet = limits(2);

p = 1;
fR = x./(1+p*(1-x));
Zs = Zvac*(Zmet/Zvac).^fR;


if nargout > 1 % provide gradient on-fly
    dfR = (p+1)./(1+p*(1-x)).^2;
    dZs = dfR .* Zvac * log(Zmet/Zvac) .* (Zmet/Zvac) .^ fR;
end
