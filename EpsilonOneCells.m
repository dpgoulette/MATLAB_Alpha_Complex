function cells1 = EpsilonOneCells(DT,VV,VC,Edges)
%

cells1 = zeros(length(Edges),3);
for i=1:size(Edges,1)
    B=intersect(VC{Edges(i,1)}, VC{Edges(i,2)});
    Poly=VV(B,:);
    edge=DT.X(Edges(i,:),:);
%             scatter3(Poly(:,1),Poly(:,2),Poly(:,3))
%             hold on
%             axis equal
%             scatter3(edge(:,1),edge(:,2),edge(:,3),'r*')
%             axis tight
%             hold off
    Poly2=bsxfun(@minus,Poly,edge(1,:));
    e2=edge(2,:)-edge(1,:);
    %figure(2)
%         scatter3(Poly2(:,1),Poly2(:,2),Poly2(:,3))
%         hold on
%         axis equal
%         scatter3([0;e2(1)],[0;e2(2)],[0;e2(3)],'r*')
%         axis tight
%         hold off
    NormalVec = [0 e2(3) -e2(2)];%cross product of e2 and [1 0 0]
    Theta=acos(e2(1)/norm(e2));%angle between e2 and [1 0 0]
    RotMat = Rodrigues(Theta,NormalVec);
    PolyRot = RotMat*(Poly2)'; %rotate the polygon. Now parallel to y-z plane
    %Just for testing
%     e2rot=RotMat*e2';
%     
%         %     figure(3)
%         scatter3(PolyRot(1,:)',PolyRot(2,:)',PolyRot(3,:)')
%         hold on
%         axis equal
%         scatter3([0;e2rot(1)],[0;e2rot(2)],[0;e2rot(3)],'r*')
%         axis tight
%         hold off
    
    PolyRot(1,:)=[];   %project onto the y-z plane. Poly is Kx2 now.
    PolyRot=PolyRot';
    if ismember(PolyRot,[0 0],'rows') %midpoint was a vertex of Poly
        fprintf('Midpoint of edge %d was a vertex of its polygon!!!\n.',i)
        M= (edge(1,:)-edge(2,:))*(1/2);
        d= norm(edge(1,:)-M);
        cells1(i,:)=[Edges(i,:),d];
        continue
    end
    %dt=DelaunayTri(PolyRot);
    MM=[PolyRot;0 0];
    originPolyID=size(MM,1);
    K=convhull(MM);
%         plot(MM(K,1),MM(K,2),'r-',MM(:,1),MM(:,2),'bo')
%         axis equal
% 
%         %axis tight

    if ~(K==originPolyID)%1) 
        %then the midpoint of e1 and e2 was in the polygon to begin with
        %so we can get the ponit directly
        M= (edge(1,:)+edge(2,:))*(1/2);%midpoint
        d= norm(edge(1,:)-M);
        cells1(i,:)=[Edges(i,:),d];
        
%             hold on
%             plot(0,0,'ro',0,0,'r+')
%             hold off
        
%             scatter3(Poly(:,1),Poly(:,2),Poly(:,3))
%             hold on
%             axis equal
%             scatter3(edge(:,1),edge(:,2),edge(:,3),'r*')
%             axis tight
%             scatter3(M(1,1),M(1,2),M(1,3),'ro')
%             scatter3(M(1,1),M(1,2),M(1,3),'r+')
%             hold off

    elseif size(MM,1)==size(K,1)-1 %then MM was convex. No verts missed.
        %don't need to do new convhull
        Oind=find(K==originPolyID,1,'first');
              
        Vbefore=MM(K(Oind-1),:);%preceding vertex in K
        Vafter =MM(K(Oind+1),:);%following vertex in K
        Nbefore = norm(Vbefore);
        Nafter = norm(Vafter);
        if Nbefore < Nafter
            V1 = Vbefore-Vafter;
            Angle = VecAngle(V1,Vbefore);
            if Angle > pi/2
                    
%                     K2=convhull(PolyRot);
%                     plot(PolyRot(K2,1),PolyRot(K2,2),'k-',PolyRot(:,1),PolyRot(:,2),'bo')
%                     axis equal
%                     hold on
%                     plot(0,0,'go')
%                     plot(Vbefore(1),Vbefore(2),'ro')
%                     plot(Vbefore(1),Vbefore(2),'r+')
%                     hold off
                
                d=norm(edge(1,:)-Poly(K(Oind-1),:));
                cells1(i,:)=[Edges(i,:),d];%The most important line!
                
%                     scatter3(Poly(:,1),Poly(:,2),Poly(:,3))
%                     hold on
%                     axis equal
%                     scatter3(edge(:,1),edge(:,2),edge(:,3),'r*')
%                     axis tight
%                     scatter3(Poly(K(Oind-1),1),Poly(K(Oind-1),2),Poly(K(Oind-1),3),'ro')
%                     scatter3(Poly(K(Oind-1),1),Poly(K(Oind-1),2),Poly(K(Oind-1),3),'r+')
%                     hold off 
            else %Angle <= pi/2
                %the case when vbefore is closer with acute angle
                
                %Angle is not obtuse so the nearest is in the
                %middle of the edge. Distance from Vbefor to the nearest point,
                % is L.
                L = Nbefore*cos(Angle);
                t=L/norm(V1);

%                 %for testing
%                     K2=convhull(PolyRot);
%                     plot(PolyRot(K2,1),PolyRot(K2,2),'k-',PolyRot(:,1),PolyRot(:,2),'bo')
%                     axis equal
%                     hold on
%                     plot(0,0,'go')
%                     N=PolyRot(K(Oind-1),:)+(PolyRot(K(Oind+1),:)-PolyRot(K(Oind-1),:))*t;
%                     plot(N(1),N(2),'ro',N(1),N(2),'r+');
%                     hold off

                %find the real point in 3d corresponding to n in 2d. Then
                %find the distance from e1 to this point.  We do this with
                %a paramaterized line in 3D
                n3=Poly(K(Oind-1),:)+(Poly(K(Oind+1),:)-Poly(K(Oind-1),:))*t;
                d=norm(edge(1,:)-n3);
                cells1(i,:)=[Edges(i,:),d];

%                 %for testing
%                     scatter3(Poly(:,1),Poly(:,2),Poly(:,3))
%                     hold on
%                     axis equal
%                     scatter3(edge(:,1),edge(:,2),edge(:,3),'r*')
%                     axis tight
%                     scatter3(n3(1),n3(2),n3(3),'ro')
%                     scatter3(n3(1),n3(2),n3(3),'r+')
%                     hold off
            end
        else %Nbefore > Nafter
            V1 = Vafter-Vbefore;
            Angle = VecAngle(V1,Vafter);
            if Angle > pi/2
                    
%                     K2=convhull(PolyRot);
%                     plot(PolyRot(K2,1),PolyRot(K2,2),'k-',PolyRot(:,1),PolyRot(:,2),'bo')
%                     axis equal
%                     hold on
%                     plot(0,0,'go')
%                     plot(Vafter(1),Vafter(2),'ro')
%                     plot(Vafter(1),Vafter(2),'r+')
%                     hold off
                
                d=norm(edge(1,:)-Poly(K(Oind+1),:));
                cells1(i,:)=[Edges(i,:),d];
                
%                     scatter3(Poly(:,1),Poly(:,2),Poly(:,3))
%                     hold on
%                     axis equal
%                     scatter3(edge(:,1),edge(:,2),edge(:,3),'r*')
%                     axis tight
%                     scatter3(Poly(K(Oind+1),1),Poly(K(Oind+1),2),Poly(K(Oind+1),3),'ro')
%                     scatter3(Poly(K(Oind+1),1),Poly(K(Oind+1),2),Poly(K(Oind+1),3),'r+')
%                     hold off  
            else %Angle <= pi/2
            %Angle is not obtuse so the nearest is in the
            %middle of the edge. Distance from Vafter to the nearest point,
            % is L.
            L = Nafter*cos(Angle);
            t=L/norm(V1);
            
%                 %for testing
%                 K2=convhull(PolyRot);
%                 plot(PolyRot(K2,1),PolyRot(K2,2),'k-',PolyRot(:,1),PolyRot(:,2),'bo')
%                 axis equal
%                 hold on
%                 plot(0,0,'go')
%                 N=PolyRot(K(Oind+1),:)+(PolyRot(K(Oind-1),:)-PolyRot(K(Oind+1),:))*t;
%                 plot(N(1),N(2),'ro',N(1),N(2),'r+');
%                 hold off
            
            %find the real point in 3d corresponding to n in 2d. Then
            %find the distance from e1 to this point.  We do this with
            %a paramaterized line in 3D
            n3=Poly(K(Oind+1),:)+(Poly(K(Oind-1),:)-Poly(K(Oind+1),:))*t;
            d=norm(edge(1,:)-n3);
            cells1(i,:)=[Edges(i,:),d];
            
            %for testing
%                 scatter3(Poly(:,1),Poly(:,2),Poly(:,3))
%                 hold on
%                 axis equal
%                 scatter3(edge(:,1),edge(:,2),edge(:,3),'r*')
%                 axis tight
%                 scatter3(n3(1),n3(2),n3(3),'ro')
%                 scatter3(n3(1),n3(2),n3(3),'r+')
%                 hold off
            
            end
        end
        
    else %then MM missed some of the poly vertices vertices.
        %find the index of the origin in K
        Oind=find(K==originPolyID,1,'first');
        
        %The indices nbring origin in K/MM.  These are the start and end of
        %a chain in Polyrot/K2
        BeginningID = K(Oind-1);%index of MM with preceding nbr of origin
        Poly=circshift(Poly,-(BeginningID-1));
        PolyRot=circshift(PolyRot,-(BeginningID-1));
        
        %EndID = K(Oind+1);%shouldn't need this
        
        Missed=size(MM,1)-size(K,1)+1;
        K2=convhull(PolyRot);
        
%             plot(PolyRot(K2,1),PolyRot(K2,2),'k-',PolyRot(:,1),PolyRot(:,2),'bo',0,0,'go')
%             axis equal
%             hold on
            %plot(0,0,'go')
            
            
        WallIDS=1:2+Missed;%indices of K that face the origin
        
        Wall=PolyRot(K2(WallIDS),:);%poly wall facing [0 0]
        
%             plot(Wall(:,1),Wall(:,2),'yo')
        
        
        WallNorms=sqrt(sum(Wall.^2,2));
        [WallMinD,WallMinID]=min(WallNorms);
        
        if WallMinID == 1
            V1=Wall(1,:)-Wall(2,:);
            Angle=VecAngle(V1,Wall(1,:));
            if Angle > pi/2
                %then wall(1) is closest
                
%                     plot(Wall(1,1),Wall(1,2),'ro')
%                     plot(Wall(1,1),Wall(1,2),'r+')
%                     hold off
                    
                d=norm(edge(1,:)-Poly(1,:));
                cells1(i,:)=[Edges(i,:),d];
                
%                     scatter3(Poly(:,1),Poly(:,2),Poly(:,3))
%                     hold on
%                     axis equal
%                     scatter3(edge(:,1),edge(:,2),edge(:,3),'r*')
%                     axis tight
%                     scatter3(Poly(1,1),Poly(1,2),Poly(1,3),'ro')
%                     scatter3(Poly(1,1),Poly(1,2),Poly(1,3),'r+')
%                     hold off
            else
                L = WallMinD*cos(Angle);
                t = L/norm(V1);
                
%                     %for testing
%                     N=PolyRot(1,:)+(PolyRot(K2(2),:)-PolyRot(1,:))*t;
%                     plot(N(1),N(2),'ro',N(1),N(2),'r+');
%                     plot(N(1),N(2),'ro',N(1),N(2),'ro');
%                     hold off
            
                %find the real point in 3d corresponding to n in 2d. Then
                %find the distance from e1 to this point.  We do this with
                %a paramaterized line in 3D
                n3=Poly(1,:)+(Poly(K2(2),:)-Poly(1,:))*t;
                d=norm(edge(1,:)-n3);
                cells1(i,:)=[Edges(i,:),d];
            
%                     %for testing
%                     scatter3(Poly(:,1),Poly(:,2),Poly(:,3))
%                     hold on
%                     axis equal
%                     scatter3(edge(:,1),edge(:,2),edge(:,3),'r*')
%                     axis tight
%                     scatter3(n3(1),n3(2),n3(3),'ro')
%                     scatter3(n3(1),n3(2),n3(3),'r+')
%                     hold off
            end
                
        elseif WallMinID == length(WallIDS)
            V1=Wall(end,:)-Wall(end-1,:);
            Angle=VecAngle(V1,Wall(end,:));
            if Angle > pi/2
                %then wall(end) is closest
                
%                     plot(Wall(end,1),Wall(end,2),'ro')
%                     plot(Wall(end,1),Wall(end,2),'r+')
%                     hold off
                    
                d=norm(edge(1,:)-Poly(K2(WallIDS(end)),:));
                cells1(i,:)=[Edges(i,:),d];
                
%                     scatter3(Poly(:,1),Poly(:,2),Poly(:,3))
%                     hold on
%                     axis equal
%                     scatter3(edge(:,1),edge(:,2),edge(:,3),'r*')
%                     axis tight
%                     scatter3(Poly(K2(WallIDS(end)),1),Poly(K2(WallIDS(end)),2),Poly(K2(WallIDS(end)),3),'ro')
%                     scatter3(Poly(K2(WallIDS(end)),1),Poly(K2(WallIDS(end)),2),Poly(K2(WallIDS(end)),3),'ro')
%                     hold off
            else
                L = WallMinD*cos(Angle);
                t = L/norm(V1);
                
                    %for testing
%                     N=PolyRot(K2(WallIDS(end)),:)+(PolyRot(K2(WallIDS(end-1)),:)-PolyRot(K2(WallIDS(end)),:))*t;
%                     plot(N(1),N(2),'ro',N(1),N(2),'r+');
%                     plot(N(1),N(2),'ro',N(1),N(2),'ro');
%                     hold off
            
                %find the real point in 3d corresponding to n in 2d. Then
                %find the distance from e1 to this point.  We do this with
                %a paramaterized line in 3D
                n3=Poly(K2(WallIDS(end)),:)+(Poly(K2(WallIDS(end-1)),:)-Poly(K2(WallIDS(end)),:))*t;
                d=norm(edge(1,:)-n3);
                cells1(i,:)=[Edges(i,:),d];
            
%                     %for testing
%                     scatter3(Poly(:,1),Poly(:,2),Poly(:,3))
%                     hold on
%                     axis equal
%                     scatter3(edge(:,1),edge(:,2),edge(:,3),'r*')
%                     axis tight
%                     scatter3(n3(1),n3(2),n3(3),'ro')
%                     scatter3(n3(1),n3(2),n3(3),'r+')
%                     hold off
            end
        else %the nearest wall vert is in the middle of the wall
            V1=Wall(WallMinID,:)-Wall(WallMinID-1,:);
            V2=Wall(WallMinID,:)-Wall(WallMinID+1,:);
            Angle1=VecAngle(V1,Wall(WallMinID,:));
            Angle2=VecAngle(V2,Wall(WallMinID,:));
            
            if Angle1 > pi/2 && Angle2 > pi/2
                %then Wall(WallMinID,:) is closest
                
%                     plot(Wall(WallMinID,1),Wall(WallMinID,2),'ro')
%                     plot(Wall(WallMinID,1),Wall(WallMinID,2),'r+')
%                     hold off
                    
                d=norm(edge(1,:)-Poly(K2(WallIDS(WallMinID)),:));
                cells1(i,:)=[Edges(i,:),d];
                
%                     scatter3(Poly(:,1),Poly(:,2),Poly(:,3))
%                     hold on
%                     axis equal
%                     scatter3(edge(:,1),edge(:,2),edge(:,3),'r*')
%                     axis tight
%                     scatter3(Poly(K2(WallIDS(WallMinID)),1),Poly(K2(WallIDS(WallMinID)),2),Poly(K2(WallIDS(WallMinID)),3),'ro')
%                     scatter3(Poly(K2(WallIDS(WallMinID)),1),Poly(K2(WallIDS(WallMinID)),2),Poly(K2(WallIDS(WallMinID)),3),'r+')
%                     hold off
            elseif Angle1 < pi/2
                L = WallMinD*cos(Angle1);
                t = L/norm(V1);
                
%                     %for testing
%                     N=PolyRot(K2(WallIDS(WallMinID)),:)+(PolyRot(K2(WallIDS(WallMinID-1)),:)-PolyRot(K2(WallIDS(WallMinID)),:))*t;
%                     plot(N(1),N(2),'ro',N(1),N(2),'r+');
%                     plot(N(1),N(2),'ro',N(1),N(2),'ro');
%                     hold off
            
                %find the real point in 3d corresponding to n in 2d. Then
                %find the distance from e1 to this point.  We do this with
                %a paramaterized line in 3D
                n3=Poly(K2(WallIDS(WallMinID)),:)+(Poly(K2(WallIDS(WallMinID-1)),:)-Poly(K2(WallIDS(WallMinID)),:))*t;
                d=norm(edge(1,:)-n3);
                cells1(i,:)=[Edges(i,:),d];
            
%                     %for testing
%                     scatter3(Poly(:,1),Poly(:,2),Poly(:,3))
%                     hold on
%                     axis equal
%                     scatter3(edge(:,1),edge(:,2),edge(:,3),'r*')
%                     axis tight
%                     scatter3(n3(1),n3(2),n3(3),'ro')
%                     scatter3(n3(1),n3(2),n3(3),'r+')
%                     hold off
            else
                L = WallMinD*cos(Angle2);
                t = L/norm(V2);
                
%                     %for testing
%                     N=PolyRot(K2(WallIDS(WallMinID)),:)+(PolyRot(K2(WallIDS(WallMinID+1)),:)-PolyRot(K2(WallIDS(WallMinID)),:))*t;
%                     plot(N(1),N(2),'ro',N(1),N(2),'r+');
%                     plot(N(1),N(2),'ro',N(1),N(2),'ro');
%                     hold off
            
                %find the real point in 3d corresponding to n in 2d. Then
                %find the distance from e1 to this point.  We do this with
                %a paramaterized line in 3D
                n3=Poly(K2(WallIDS(WallMinID)),:)+(Poly(K2(WallIDS(WallMinID+1)),:)-Poly(K2(WallIDS(WallMinID)),:))*t;
                d=norm(edge(1,:)-n3);
                cells1(i,:)=[Edges(i,:),d];

            end
        end
    end
end
end

function RotMat = Rodrigues(Theta,u)
% Rodrigues takes a vector axis of rotation and an angle of rotation. The
%       3-d rotation matrix is returned.  This is Rodrigues' rotation
%       formula.

u = u./norm(u,2);
S = [    0  -u(3) u(2);
    u(3)   0   -u(1);
    -u(2) u(1)   0  ];
RotMat = eye(3) + sin(Theta)*S + (1-cos(Theta))*S^2;

end

function Theta = VecAngle(v1,v2)
%  find the angle between two row vectors

Theta = acos((v1*v2')/(norm(v1)*norm(v2)));
end
