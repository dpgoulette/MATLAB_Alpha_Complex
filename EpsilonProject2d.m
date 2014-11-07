DT=DelaunayTri(Flatdata);
[VV VC]=voronoiDiagram(DT);

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

% %Now we can use this BadData vector to remove any Triangles and
% %Edges that contain any of these bad points.  Then we will only be
% %considering cells with GoodPoints for inclusion to the epsilon complex.

%get ALL tris and edges (including bad ones)
DTtris=DT.Triangulation;
DTedges=edges(DT);


%now remove bad triangles.  They are "bad" if they contain
%a BadData point.

[R,~]=find(ismember(DTtris,BadData));
R=unique(R);
cells2=DTtris;
cells2(R,:)=[];
%find the circumcenter radius for the triangles.  This is epsilon.
[~,RCC]=circumcenters(DT);
RCC(R)=[];%delete the CircCenters of bad triangles

cells2=[cells2, RCC];%2-cells are done.



[R,~]=find(ismember(DTedges,BadData));
R=unique(R);
goodEdges=DTedges;
goodEdges(R,:)=[];
%find the epsilon for each edge

cells1 = EpsilonOneCells2d(DT.X,goodEdges,VV,VC);

clear a b R ID BadVV goodEdges