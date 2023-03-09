function Zrho = ConstructMaterialMatrix(Mesh, BF, Zs, triangles)
%% Construct full material matrix Zrho for MoM eq. or a element material matrix
%
% Inputs:
%   Mesh        - MATLAB structure of geometry (from AToM, see [1])
%   BF          - MATLAB structure of basis functions (from AToM, see [1])
%   Zs          - surface resistivities of triangles 
%   triangles   - unambiguous labels of triangles for which is Zrho
%                 computed. If triangles=1:T the full material matrix is
%                 evaluated.  
%
% Outputs:
%   Zrho        - material matrix [N x N]
% 
% 2023, Jonas Tucek, CTU in Prague, jonas.tucek@fel.cvut.cz

%% prealocation
nodes = Mesh.nodes;
nTriangles = Mesh.nTriangles;
edgeLengths = Mesh.triangleEdgeLengths;
trAreas = Mesh.triangleAreas;
trCentroids = Mesh.triangleCentroids;
trNodes = Mesh.connectivityList;
nUnknowns = BF.nUnknowns;
BFdata = BF.data;

if nargin < 3 || isempty(Zs)
   Zs = ones(nTriangles, 1);
end

% material matrix
Zrho  = zeros(nUnknowns);

for iTriangle = triangles
   % edges at plus triangle
   edgesAtTriaPlus = BFdata(BFdata(:, 2) == iTriangle, 7);
   % edges at minus triangle
   edgesAtTriaMinus = BFdata(BFdata(:, 4) == iTriangle, 7);

   % aggregate all valid edges
   edgesAtTria = [edgesAtTriaPlus; edgesAtTriaMinus];
   nEdgesAtTria = length(edgesAtTria);

   % basis functions' orientation
   xi = [ones(length(edgesAtTriaPlus), 1); ...
      -1*ones(length(edgesAtTriaMinus), 1)];

   % length of all relevant edges
   edgeLengthsAtTria = edgeLengths(BFdata(edgesAtTria, 3));

   % geometry constants for this triangle
   Rt = trCentroids(iTriangle, :);
   Pt = nodes(trNodes(iTriangle, :), :);
   RPcoef = 3/4 * dot(Rt, Rt) + 1/12 * sum(dot(Pt, Pt, 2));
   
   % loop over all edges
   edgeCounter = 1;
   for iEdge1 = 1:nEdgesAtTria
      % free vertex to first edge
      % depend on plus (column 1) or minus (column 2) triangle:
      % column = 3 - 2*xi(iEdge)
      P1n = nodes(BFdata(edgesAtTria(iEdge1), 3 - 2*xi(iEdge1)), :);
      for iEdge2 = 1:nEdgesAtTria
         % free vertex to second edge
         P2n = nodes(BFdata(edgesAtTria(iEdge2), 3 - 2*xi(iEdge2)), :);

         % contribution to lossy matrix
         Litem = 1/4 *  xi(iEdge1) * xi(iEdge2) / trAreas(iTriangle) * ...
            edgeLengthsAtTria(iEdge1) * edgeLengthsAtTria(iEdge2) * ...
            (RPcoef - dot(Rt, P1n + P2n) + dot(P1n, P2n));

         % save this entry to the final L matrix, apply triangle resitivity
         Zrho(edgesAtTria(iEdge1), edgesAtTria(iEdge2)) = ...
            Zrho(edgesAtTria(iEdge1), edgesAtTria(iEdge2)) + ...
            Zs(iTriangle) * Litem;
         edgeCounter = edgeCounter + 1;
      end
   end
end
end

