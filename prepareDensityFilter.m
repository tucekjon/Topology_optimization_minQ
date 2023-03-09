function H = prepareDensityFilter(Mesh, settings)
%% Prepare a density filter with radius settings.rmin represented by matrix H
% Filtered density field is obtained as xPhys = (H * x) ./ sum(H,2)
%
% Inputs:
%   Mesh        - MATLAB structure of geometry (from AToM, see [1])
%   settings    - MATLAB structure containing optimization settings.
%                 Namely, radius of density filter Rmin, protected
%                 triangles and passive triangles labels, are required.
%
% Outputs:
%   H      - matrix representing a convolution operation of density filtering (see the paper)
%
% 2023, Jonas Tucek, CTU in Prague, jonas.tucek@fel.cvut.cz


H = zeros(Mesh.nTriangles,Mesh.nTriangles);

for iT = 1:Mesh.nTriangles
    rVec = Mesh.triangleCentroids(iT,:) - Mesh.triangleCentroids; % Radius vectors of iT center to all other triangles
    dist = sqrt(sum(rVec.^2,2));     % Distance iT triangles to all other triangles
    indDist = dist <= settings.rmin; % Labels of neighboring triangles

    H(iT,indDist) = settings.rmin - dist(indDist);
end

H = sparse(H);

% Do not perform density filtering on protected triangles
nProtectedTRs = length(settings.protTRs);
H(settings.protTRs,:) = 0;
H(:,settings.protTRs) = 0;
H(settings.protTRs,settings.protTRs) = eye(nProtectedTRs);

% Do not perform density filtering on passive triangles
nPassiveTRs = length(settings.passiveTRs);
H(settings.passiveTRs,:) = 0;
H(:,settings.passiveTRs) = 0;
H(settings.passiveTRs,settings.passiveTRs) = eye(nPassiveTRs);
end