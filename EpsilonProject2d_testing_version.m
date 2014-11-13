DT=DelaunayTri(FlatDataExample);
[VV VC]=voronoiDiagram(DT);

% In our work we choose to remove data that is on the boundary of the data
% space.  We do this by removing any data point that has a voronoi cell
% vertex which we deem to be a "bad Voronoi vertex."  We also remove any
% Infinity is bad because any voronoi cell atteched to it has infinite
% area/volume.  Also any voronoi vertex outside of the data space is bad
% because they are often artificially elongated, thus the voronoi cell is
% less reliable.

% Note, any finite voronoi vertex that is outside of the dataspace
% will not be inside any Delaunay tetra.

% find the indices of the bad Voronoi vertices.  pointLocation returns a nan
% if the voronoi vertex is NOT interior to any triangle/tetra.  We put a NaN
% in the first slot because the point at infinity is listed first in the
% Voronoi diagram.
ID = [NaN; pointLocation(DT,VV(2:length(VV),:))];
BadVV = find(isnan(ID));%creates index vector of bad vor verts

% Now we find the bad data.  "Bad data" is on the boundary of the data space
% so it has ind the vertices in the data that have a voronoi cell which
% includes a bad voronoi vertex.  Create an index of "BadDataID."
BadDataID = zeros(length(DT.X),1);
a=1;
for b=1:length(VC)
    if ~isempty(intersect(BadVV,VC{b}))
        BadDataID(a)=b;
        a=a+1;
    end
end
BadDataID(BadDataID==0)=[];

% Now we can use this BadData vector to remove any Triangles and
% edges that contain any of these bad points.  Then we will only be
% considering cells with GoodPoints for inclusion to the epsilon complex.

% Get ALL tris and edges (including bad ones).  These are the 1 and 2
% cells. 
DTtris=DT.Triangulation;
DTedges=edges(DT);

% Now remove bad triangles.  They are "bad" if they contain

[R,~] = find(ismember(DTtris, BadDataID));
R = unique(R);
cells2 = DTtris;
cells2(R,:) = [];

% Find the circumcenter radius for the triangles.  This is the value of
% alpha for all of teh 2 cells.
[~,RCC] = circumcenters(DT);
% Now delete the radii from bad triangles
RCC(R)=[];
% Save the 2 cells along with their alpha value in the third column.
cells2=[cells2, RCC];

% Remove the bad one cells (edges).
[R,~]=find(ismember(DTedges,BadDataID));
R=unique(R);
goodEdges=DTedges;
goodEdges(R,:)=[];
% Find the epsilon for each edge
cells1 = EpsilonOneCells2d(DT,goodEdges,VV,VC);

clear a b R ID BadVV goodEdges