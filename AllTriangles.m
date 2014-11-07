function Triangles = AllTriangles(DT)
%

%MAKE THE MATRIX OF UNIQUE TRIANGLES
Triangles=zeros(4*length(DT.Triangulation),3);
for a=1:length(DT.Triangulation)
    Triangles((a-1)*4+1:(a-1)*4+1+3,:)=nchoosek(DT.Triangulation(a,:),3);
    fprintf('%d finished.\n',a)
end
end