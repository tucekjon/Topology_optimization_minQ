function history = topOptInMoM(OP, settings)
%% Perform the denstity topology optimization in Method of Moments Models 
%   Antenna Q-factor is minimized by distributing a conductive material
%   Refer to the paper [1] for more details.
% 
% Inputs:
%   OP           - MATLAB structure containing all variables and fields
%                  fully describing the optimization region and all
%                  necessary MoM matrices. Majority of the fields are
%                  evaluated by AToM, see [1]. Mandatory fields are:
%                  OP.Mesh (discretization), OP.BF (basis functions),
%                  OP.ports (position of discrete ports), OP.V
%                  (excitation vector), OP.Z0 (vacuum impedance matrix),                  
%                  etc. See START.m for details.
%   settings     - MATLAB structure containing optimization settings.
%                  See START.m for details.
% 
% Outputs:
%   history      - MATLAB structure containing the history of the optimization process.
% 
% The code is initiated from START.m.
% 
% 2023, Jonas Tucek, CTU in Prague, jonas.tucek@fel.cvut.cz

%% Prepare Density filter
H = prepareDensityFilter(OP.Mesh,settings);

%% Initialize iteration
x = settings.Sf*ones(OP.Mesh.nTriangles,1); % Initial design field distribution
x(settings.protTRs) = 1; % Protected triangles are metal
x(settings.passiveTRs) = 0; % Passive triangles are vacuum

% Indices of design variables (triangles)
settings.designTRs = setdiff(1:OP.Mesh.nTriangles,[settings.passiveTRs settings.protTRs]);

% Total design domain surface
S0 = sum(OP.Mesh.triangleAreas(settings.designTRs));

% Initial plot, which is to be updated
handle = plotDesign(OP.Mesh, OP.BF, x, OP.port);

change = 1;
iter = 0;
loopbeta = 0;
nProtectedTRs = length(settings.protTRs);
nPassiveTRs = length(settings.passiveTRs);

%% Setup varibles for History struct
history.x = x;                  % Design variable
history.fval = [];              % Q-factors values
history.change = [];            % change between consecutive iterations
history.beta = settings.beta;   % sharpness of the projection filter
history.constrVal = [];         % value of the are fraction constraint

%% MMA setup
move = 0.25; % External moving limits
m     = 3;                             % The number of general constraints.
n     = OP.Mesh.nTriangles...
    - nProtectedTRs - nPassiveTRs; % The number of design variables x_j.
xmin  = zeros(n,1);                    % Column vector with the lower bounds for the variables x_j.
xmax  = ones(n,1);                     % Column vector with the upper bounds for the variables x_j.
xold1 = x(settings.designTRs);         % xval, one iteration ago (provided that iter>1).
xold2 = x(settings.designTRs);         % xval, two iterations ago (provided that iter>2).
low   = zeros(n,1);                    % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp   = ones(n,1);                     % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
a0    = 1;                             % The constants a_0 in the term a_0*z.
a     = [ones(m-1,1);0];               % Column vector with the constants a_i in the terms a_i*z.
c     = 1e3*ones(m,1);                 % Column vector with the constants c_i in the terms c_i*y_i.
d     = zeros(m,1);                    % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.

%% Start the topology optimization task

% Initial projection filter (with beta = 1 has zero influence)
[xPhys, dx] = projectionFilter(x, settings.beta, settings.etaVec);

tStart=tic;
while change > settings.change && iter < settings.maxIter
    iter = iter + 1;         % Iteration
    loopbeta = loopbeta + 1; % Number of iterations for projection filter update

    % Evaluate objectives and Adjoint sensitivity analysis
    [f, df] = settings.fitness(OP, settings, xPhys);

    % Use chain-rule to determine sensitivities with respect to the design variable x
    df = H * (df .* dx ./sum(H,2)); % Senstivities of objectives

    dS = OP.Mesh.triangleAreas;      % Areas of each triangle
    dS = H * (dS .* dx ./sum(H,2)); % Sensitivities of Area fraction constraint
    %% MMA update of the design variable
    % MMA is set according to the optimization task (see [1])
    xval = x(settings.designTRs); % Only design triangles are updated
    % Objective value and sensitivities
    f0val = 0;
    df0dx = 0;
    % Constraint values and sensitivities
    fval = [f(:);...
        (xPhys(settings.designTRs).' * OP.Mesh.triangleAreas(settings.designTRs))/(settings.Sf*S0) - 1];
    dfdx = [df(settings.designTRs,:).';...
        dS(settings.designTRs).'/(settings.Sf*S0)];

    % Perform MMA update
    [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
        mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
    
    % Retrieve updated design variable
    xnew = ones(OP.Mesh.nTriangles,1);
    xnew(settings.passiveTRs) = 0;
    xnew(settings.designTRs) = xmma;
    
    % External moving limit update
    xmin  = max(xnew(settings.designTRs)-move,zeros(n,1));
    xmax  = min(xnew(settings.designTRs)+move,ones(n,1));

    % Perform filtering
    xTilde = (H * xnew) ./sum(H,2); % Density filter
    [xPhys, dx] = projectionFilter(xTilde, settings.beta, settings.etaVec); % Projection filter

    xold2    = xold1;                 % Saved for next MMA update
    xold1    = x(settings.designTRs); % Saved for next MMA update
    change = max(abs(xnew-x));        % Change between consecutive runs
    x = xnew;                         % New design variable

    %% Print and plot
    TrianglesColour = ([255 255 255] - xPhys*([255 255 255] - [0 0 0]))/255;
    handle.FaceVertexCData = TrianglesColour;
    pause(0.01);
    fprintf('iter.:%5i AreaFrac.:%7.3f change:%7.3f Qe:%1.3f, Qm:%1.3f\n',iter, ...
        (xPhys(settings.designTRs).' * OP.Mesh.triangleAreas(settings.designTRs))/(S0), change, f(1), f(2));

    %% Save History
    history.constrVal = [history.constrVal (xPhys(settings.designTRs).' * OP.Mesh.triangleAreas(settings.designTRs))/S0];
    history.beta = [history.beta settings.beta];
    history.x = [history.x x];
    history.fval = [history.fval f(:)];
    history.change = [history.change change];

    %% Update sharpness of the projection filter
    if settings.beta < settings.betaMax && (loopbeta >= settings.betaIter || change <= settings.change)
        settings.beta = 2 * settings.beta; % Double sharpness value
        change = 1;
        loopbeta = 0;
        fprintf('Sharpness of the projection filter, beta, increased to %g.\n',settings.beta);
    end
end
t0 = toc(tStart);
fprintf(['\nTopology optimization found the design with Q-factor=%1.3f in %3i iterations and %1.2f seconds occupying  %1.1f%%' ...
    ' of the design domain\n'],max(f,[],2), iter, t0, history.constrVal(end)*100); 

history.t0 = t0;
history.H = H;
