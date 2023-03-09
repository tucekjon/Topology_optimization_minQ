function [xOut, dxOut] = projectionFilter(xIn, beta, eta)
% Smooth approximation of Heaviside step function
% - Performs smooth thresholding with sharpness beta centered
%   at eta.
% 
% Inputs:
%   xIn             - variable to be projected by the filter
%   beta            - sharpness of the filter
%   eta             - level of the filter
% 
% Outputs:
%   xOut      - projected variable
%   dxOut     - derivative of the projection filter
% 
% 2023, Jonas Tucek, CTU in Prague, jonas.tucek@fel.cvut.cz

% Projection filter
xOut = (tanh(beta*eta)+tanh(beta*(xIn-eta)))./(tanh(beta*eta)+tanh(beta*(1-eta)));

if nargout > 1 % Provide derivative on-fly
    dxOut = (1-tanh(beta*(xIn-eta)).^2)*beta./(tanh(beta*eta)+tanh(beta*(1-eta))); 
end

