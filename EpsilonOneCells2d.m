function cells1 = EpsilonOneCells2d(DT,Edges,VV,VC)
%

Data = DT.X;
cells1 = [Edges, zeros(length(Edges),1)];

for i=1:size(Edges,1)
    B=intersect(VC{Edges(i,1)}, VC{Edges(i,2)});
   
    if size(B,2) ~= 2
        fprintf('ERROR! This was unexpected.\n')
        error('Neighboring data-points in 2D must have two Voronoi vertices in common.')
    else
        % Pull one of the data points from the current Delaunay edge.  (Both
        % endpoints are equidistant to the voronoi wall so it doesn't
        % matter which one.)
        x=Data(Edges(i,1),:);
        % Pull one endpoint of the voronoi edge they share
        e1=VV(B(1),:);
        % Pull the other endpoint of the voronoi edge they share
        e2=VV(B(2),:);
        
        v1=e2-e1;
        v2=x-e1;
        Theta1 = acos((v1*v2')/(norm(v1)*norm(v2)));
        
        w1=e1-e2;
        w2=x-e2;
        Theta2 = acos((w1*w2')/(norm(w1)*norm(w2)));
        
        if Theta1 > pi/2 %then e1 is closest point
            cells1(i,3) = norm(v2);
        elseif Theta2 > pi/2 %then e2 is closest point
            cells1(i,3) = norm(w2);
        else
            % The distance to the midpoint is closest.  So half the edge
            % length is alpha.  This means that the Delaunay edge between
            % the data points intersects the voronoi wall they share.
            cells1(i,3) = (1/2)*norm(Data(Edges(i,1),:)-Data(Edges(i,2),:));
        end
    end
end % for loop
end % function
