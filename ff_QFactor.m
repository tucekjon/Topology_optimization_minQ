function [Q, dQ] = ff_QFactor(OP, settings, x)
%% Evalute Q-factor and compute sensitivities by Adjoint Sensitivity analysis
%       This function evaluates Q-factors, Q_e and Q_m, without any
%       penalty on the self-resonance
% 
% Inputs:
%   OP              - MATLAB structure containing all required information
%                     about the structure. Matrices Z0, Xe and Xm, BF2T are required. The
%                     excitation vector V must also be defined.
%   settings        - MATLAB structure containing optimization settings.
%                     Namely, boundary resistivities, normalization constant and design triangles labels, is necessary.
%   x               - design variable [0,1]^Tx1
% 
% Outputs:
%   Q      - interpolated surface resistivities
%   dQ     - derivative of the interpolation function  
% 
% 2023, Jonas Tucek, CTU in Prague, jonas.tucek@fel.cvut.cz

%% MoM Analysis
[Zs, dZs] = settings.interFun(x, settings.resistivityLimits); % Map design variable to resistivities
Zrho = ConstructMaterialMatrix(OP.Mesh, OP.BF, Zs, 1:OP.Mesh.nTriangles); % Construct the material matrix Zrho

Z = OP.Z0 + Zrho; % Full impedance matrix
I = Z \ OP.V;     % Evaluate current

%  Evaluate Q-factor value
R0 = real(OP.Z0);
Prad =  1/2 * I' * R0 * I;

We = 1/4 * I' * OP.Xe * I;
Wm = 1/4 * I' * OP.Xm * I;

Qe = real(2 * We/Prad); % Electric Q-factor
Qm = real(2 * Wm/Prad); % Magnetic Q-factor
Q = [Qe Qm]./settings.normalization; 

if nargout > 1 % Adjoint sensitivity analysis provided on-fly
adjR_e =  - 1/(2*Prad) * I'*(OP.Xe - Qe*R0); % Derivative of Qe to I
adjR_m =  - 1/(2*Prad) * I'*(OP.Xm - Qm*R0); % Derivative of Qm to I

lam = Z \ ([adjR_e.' adjR_m.']); % Solve two adjoint equations for each Q-factor
lam_e = lam(:,1);
lam_m = lam(:,2);

dfQe = zeros(OP.Mesh.nTriangles, 1);
dfQm = zeros(OP.Mesh.nTriangles, 1);
    
% Evaluate sensitivities
for iT = settings.designTRs    % Sweep through design triangles
    BFs = find(OP.BF2T(iT,:)); % find corresponding basis functions of the triangle
    Ie = I(BFs);               % current flowing through corresponding basis functions
    lame_e = lam_e(BFs);       
    lame_m = lam_m(BFs);
    dZrhoEl = ConstructMaterialMatrix(OP.Mesh,OP.BF, dZs, iT); % Compute element matrix
    % Final sensitivities
    dfQe(iT,:) = 2*real(lame_e.' * dZrhoEl(BFs,BFs) * Ie)./settings.normalization;
    dfQm(iT,:) = 2*real(lame_m.' * dZrhoEl(BFs,BFs) * Ie)./settings.normalization;
end

dQ = [dfQe dfQm];
end


