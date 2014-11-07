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

%find bad voronoi vertices.  Infinity is bad.  And any voronoi vertex
%outside of the data space.  Thus it is not contained in any delaunay
%tetra.
ID = [NaN; pointLocation(DT,VV(2:length(VV),:))];
BadVV = find(isnan(ID));

%find the vertices in the data that have a voronoi cell which includes a
%bad voronoi vertex
BadData = zeros(length(DT.X),1);
a=1;
for b=1:length(VC)
    if ~isempty(intersect(BadVV,VC{b}))
        BadData(a)=b;
        a=a+1;
    end
end
BadData(BadData==0)=[];

% %Now we can use this BadData vector to remove any Tetras, Triangles and
% %Edges that contain any of these bad points.  Then we will only be
% %considering cells with GoodPoints for inclusion to the epsilon complex.

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

[R,~]=find(ismember(DTtets,BadData));
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

[R,~]=find(ismember(DTtriangles,BadData));
R=unique(R);
goodTris=DTtriangles;
goodTris(R,:)=[];
cells2=EpsilonTwoCells(DT,VV,VC,goodTris);

[R,~]=find(ismember(DTedges,BadData));
R=unique(R);
goodEdges=DTedges;
goodEdges(R,:)=[];

%find the distance to their nearest common point 
%make the cells1 Jx3 matrix with good edges and min dist





