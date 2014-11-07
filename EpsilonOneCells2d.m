function cells1 = EpsilonOneCells2d(X,Edges,VV,VC)
%


cells1 = [Edges, zeros(length(Edges),1)];

for i=1:size(Edges,1)
    B=intersect(VC{Edges(i,1)}, VC{Edges(i,2)});
    
%     P=X(Edges(i,:),:);
%     V=VV(B,:);
%     plot(P(:,1),P(:,2),'k*')
%     hold on
%     axis equal
%     plot(V(:,1),V(:,2),'ko')
    
    
    if size(B,2) < 2
        fprintf('Tell David to fix this!\n')
        error('Neighbor data-points have less than one Vor Vert in common')
    else
        x=X(Edges(i,1),:);%one of the data points from the edge.
        e1=VV(B(1),:);%one endpoint of the voronoi edge they share
        e2=VV(B(2),:);%the other endpoint
        
        v1=e2-e1;
        v2=x-e1;
        Theta1 = acos((v1*v2')/(norm(v1)*norm(v2)));
        
        w1=e1-e2;
        w2=x-e2;
        Theta2 = acos((w1*w2')/(norm(w1)*norm(w2)));
        
        if Theta1 > pi/2 %then e1 is closest point
            cells1(i,3) = norm(v2);
            
%             plot(e1(1),e1(2),'rx')
            
        elseif Theta2 > pi/2 %then e2 is closest point
            cells1(i,3) = norm(w2);
            
%             plot(e2(1),e2(2),'rx')
            
        else
            %The distance to the midpoint is closest.  So half the edge
            %length is epsilon.
            cells1(i,3) = (1/2)*norm(X(Edges(i,1),:)-X(Edges(i,2),:));
            
%             M=(1/2).*(x+X(Edges(i,2),:));
%             plot(M(1),M(2),'rx')
            
        end
    end
%     hold off
end
end
