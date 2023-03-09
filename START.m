%%  Script performing density topology optimization for the Q-factor minimization 
% 
%  The method performs a local gradient-driven step evaluated by Adjoint
%  sensitivity analysis until it converges to a local minimum. The solution
%  space is regularized by the filtering technique (density and projection
%  filters) to penalize unwanted behavior. See [1] for a more detailed
%  description.
% 
%  The method uses matrices stemming from method-of-moments solution to field 
%  integral equation with Rao-Wilton-Glisson basis functions. 
%  In this script, the precalculated matrices result from 
%  electric field integral equation. Triangular mesh grid elements
%  are used as optimization degrees of freedom.
% 
%  Q-factor minimization is performed on a design domain discretized into
%  960 triangles at the electrical size ka = 0.8, see first example in [1].
% 
% [1] Tucek, J.,Capek, M., Jelinek, L., Sigmund, O.: Q-factor Minimization via Topology Optimization, 
%     pp. 1-13, 2023.

clc; 
clear; 
close all;

% Precalculated matrices are loaded (they can be changed in AToM [1]):
load('operators_plate_TopOpt.mat');

% Auxiliarly matrices
OP.Xm    = Xm;      % stored magnetic energy matrix
OP.Xe    = Xe;      % stored electric energy matrix
OP.Z0    = Z0;      % vacuum part of the impedance matrix 
OP.V     = V;       % excitation vector
OP.Mesh  = Mesh;    % MATLAB structure of geometry (from AToM)
OP.BF    = BF;      % MATLAB structure of basis functions (from AToM)   
OP.BF2T  = BF2T;    % Connectivity matrix between triangles and basis functions

%% Fitness function settings
OP.port = port;  % label of the basis function with delta gap excitation

% Q-factor fitness function 
topOptSettings.fitness = @ff_QFactor;

%% Topology optimization settings:
topOptSettings.Sf = 0.35;                           % Area fraction constraint
topOptSettings.normalization = QlbTM;               % Fitness function is normalized to fundamental bound (See FunBo [1])
topOptSettings.interFun = @interFun;                % Define interpolation function
topOptSettings.resistivityLimits = [1e5 1];         % Boundary resistivity values for vacuum (rho=0) and metal (rho=1)
topOptSettings.rmin = 0.15*OP.Mesh.normDistanceA;   % Density filter radius

% Projection filter is used in the continuation scheme
topOptSettings.beta = 1;       % Initial sharpness of the projection filter H
topOptSettings.etaVec  = 0.5;  % Level of the projection filter H
topOptSettings.betaIter = 100; % Projection filter is updated every 100 iterations
topOptSettings.betaMax = 32;   % Maximal allowed sharpness of the projection filter H

topOptSettings.change = 0.01;  % Maximal change in rho between two consecutive runs
topOptSettings.maxIter = 600;   % Maximal number of iterations

% Region near the feeder is fixed
topOptSettings.protTRs = BF.data(OP.port,[2 4]); % Fixed material on the delta gap 
topOptSettings.passiveTRs = find(abs(Mesh.triangleCentroids(:,1)) <= 1 & ...
                                 Mesh.triangleCentroids(:,2) > 0 & ...
                                 Mesh.triangleCentroids(:,2) < 4).';

% Start the optimization
history = topOptInMoM(OP, topOptSettings);

%% Postproccessing
% Convergence history
iter = 1:length(history.fval);
figure;
plot(iter, history.fval, 'x-', 'LineWidth', 1); hold on;
grid on;
xlabel('iteration');
ylabel('Q (-)');
ylim([0.5 3])
% add vertical dashed line when sharpness beta was updated
betaVec = unique(history.beta);
for i=1:length(betaVec)
    bb = find(betaVec(i) == history.beta);
    txt = ['it = ' num2str(bb(end)) ', beta = ' num2str(betaVec(i))];
    xline(bb(end),'--',{txt},'LabelHorizontalAlignment', 'left'); 
end
hold off;

legend('$Q_\mathrm{e}/Q_\mathrm{lb}^{\mathrm{TM}}$',...
       '$Q_\mathrm{m}/Q_\mathrm{lb}^{\mathrm{TM}}$',...
       '$\beta$ update','Interpreter','latex','FontSize',12)

% Plot final optimized design
iter = iter(end); % last iteration
betaPostProcess = history.beta(iter);
x = history.x(:,iter); % design variable 
xTilde = (history.H * x)./ sum(history.H,2);
xPhys = projectionFilter(xTilde, betaPostProcess, topOptSettings.etaVec); % physical (filtered) variable

plotDesign(OP.Mesh,OP.BF,xPhys,OP.port);
QGray = max(history.fval(:,end));
title(['Q-factor = ' num2str(QGray,3)])

% Hard thresholding of the gray design
xPhysBW = projectionFilter(xTilde, 150000, 0.5);
QBW = topOptSettings.fitness(OP, topOptSettings, xPhysBW);

handle = plotDesign(OP.Mesh,OP.BF,xPhysBW,OP.port); % plot structure

% Analyze BW design with physically accurate resistivites of vacuum and PEC
topOptSettings.resistivityLimits = [1e8 0.01];
Q = topOptSettings.fitness(OP, topOptSettings, xPhysBW);

title(['Q-factor = ' num2str(max(Q),3)])
