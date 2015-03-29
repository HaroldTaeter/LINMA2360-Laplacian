%
% Laplacian: 
% generate SDD matrix from electric grid
%

m=4;
file='data.txt';
myFile = fopen(file,'w+');
fprintf(myFile,'Number of nodes %d \n', m*m);
fprintf(myFile,'Number of edges %d \n', 2*m*(m-1));

A=zeros(m*m,m*m);
r=1;
D=zeros(m*m);
D(1:m)=[2, 3*r*ones(m-2), 2];
for i=1:m-2
    D(m+m*(i-1):2*m+m*(i-1))=[3, 4*r*ones(m-2), 3];
end
D(end-m:end)=[2, 3*r*ones(m-2), 2];

A=A+diag(10*ones(1,m),0)+diag(2*ones(1,m-1),1)...
    +diag(2*ones(1,m-1),-1)...
    +diag(1*ones(1,m-2),2)+diag(-1*ones(1,m-2),-2);
for i=1:m
    for j=1:m
        fprintf(myFile,' %d ', A(i,j));
    end
    fprintf(myFile,'\n');
end







