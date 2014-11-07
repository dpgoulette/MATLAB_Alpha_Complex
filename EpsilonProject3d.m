DT=DelaunayTri(flatdata);
[VV VC]=voronoiDiagram(DT);

%DTtets=DT.Triangulation;
%free=freeBoundary(DT); 

% for a=1:length(free)
%     free(a,:)=sort(free(a,:));
% end
% A=ismember(cells2,free,'rows');
% B=find(A);

% %%%%%%SAVED

% %checking the behavior of vertices that reference infinity
% test=EpsilonEdges(DT,VV,VC,cells2);
% test(282:end,:)=[];
% find(test(:,4)~=0)
% test(:,4)=[];

% %I think I need to remove any data point that has a voronoi cell that goes
% %to infinity.  Then remove any tetra that used one of these points.

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

%Now we can use this BadData vector to remove any Tetras, Triangles and
%Edges that contain any of these bad points.  Then we will only be
%considering cells with GoodPoints for inclusion to the epsilon complex.

%get ALL tets and edges (including bad ones)
DTtets=DT.Triangulation;
DTedges=edges(DT);

%get ALL triangles in DT (including bad ones)
DTtriangles=AllTriangles(DT);
for a=1:length(DTtriangles)
    DTtriangles(a,:)=sort(DTtriangles(a,:));
end
DTtriangles=unique(DTtriangles,'rows');

%now remove bad tets, edges and triangles.  They are "bad" if they contain
%a BadData point.

[R,~]=find(ismember(DTtets,BadDataID));
R=unique(R);
cells3=DTtets;
cells3(R,:)=[];

[~,RCC]=circumcenters(DT);
RCC(R)=[];%delete the CircCenters of bad tetras

cells3=[cells3, RCC];

%Now cells3 is an Nx5 matrix with the good tets and the distance to their
%nearest common vor vertex.

%find the distance to their nearest common point with EpsilonEdges
%which returns the cells2 Kx4 matrix with good tris and min dist

[R,~]=find(ismember(DTtriangles,BadDataID));
R=unique(R);
goodTris=DTtriangles;
goodTris(R,:)=[];
cells2=EpsilonTwoCells(DT,VV,VC,goodTris);

[R,~]=find(ismember(DTedges,BadDataID));
R=unique(R);
goodEdges=DTedges;
goodEdges(R,:)=[];

%find the distance to their nearest common point 
%make the cells1 Jx3 matrix with good edges and min dist





