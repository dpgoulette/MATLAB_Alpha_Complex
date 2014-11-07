function Cells2 = EpsilonTwoCells(DT,VV,VC,Tri)
%
%
%

%v is a vertex from the triangle. x and y are the vertices of the edge that
%the three data points in the triangle share.  This is the equation to find
%the parameter for the minimum distance to line that contains the edge.
%minDist =@(v,x,y) -(x-v)*(y-x)'/norm(y-x)^2;

Cells2 = zeros(length(Tri),4);
for i=1:size(Tri,1)
    B=intersect(VC{Tri(i,1)}, intersect( VC{Tri(i,2)}, VC{Tri(i,3)}));
    
    if size(B,2) == 2
        %parameter value for the min distance to the line
        %t = minDist(DT.X(Tri(i,1),:),V(B(1),:),V(B(2),:));
        
        t = -(VV(B(1),:)-DT.X(Tri(i,1),:))*(VV(B(2),:)-VV(B(1),:))'/norm(VV(B(2),:)-VV(B(1),:))^2;
        if t <= 0
            minDist = norm(DT.X(Tri(i,1),:)-VV(B(1),:));
        elseif t >= 1
            minDist = norm(DT.X(Tri(i,1),:)-VV(B(2),:));
        else
            p = VV(B(1),:)+(VV(B(2),:)-VV(B(1),:))*t;
            minDist = norm(p-DT.X(Tri(i,1),:));
        end
        Cells2(i,:) = [Tri(i,:), minDist];
        %Cells2(i,:)= [Tri(i,:),1];
        
%     elseif size(B,2) == 1
%         fprintf('Warning: a triplet of Delaunay datapoints have only one vertex in common!\n')
%         Cells2(i,:) = [Tri(i,:), -1];
%     elseif size(B,2) == 3
%         fprintf('Warning: a triplet of Delaunay datapoints has more than 2 vertices in common!\n')
%         if ~(B==1)
%             %Cells2(i,:) = [Tri(i,:), -3];
%             continue
%         else
%             Cells2(b,1:3)= B;
%             b=b+1;
%         end
    else
        fprintf('Warning: a triplet of Delaunay vertices is pathological!\n')
        Cells2(i,:) = [Tri(i,:), NaN];% NaN marks the error
    end   
    
    fprintf('%d of %d completed\n',i,size(Tri,1))
end
end
