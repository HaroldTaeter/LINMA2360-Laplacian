%
% Laplacian: 
% generate SDD matrix
%

m=5;
file='dataL.txt';
myFile = fopen(file,'w+');
fprintf(myFile,'Number of nodes %d \n', m);
fprintf(myFile,'Number of edges %d \n', (m-1)+(m-2));

A=zeros(m,m);

A=A+diag(10*ones(1,m),0)+diag(-2*ones(1,m-1),1)...
    +diag(-2*ones(1,m-1),-1)...
    +diag(-1*ones(1,m-2),2)+diag(-1*ones(1,m-2),-2);
for i=1:m
    for j=1:m
        fprintf(myFile,' %d ', A(i,j));
    end
    fprintf(myFile,'\n');
end







