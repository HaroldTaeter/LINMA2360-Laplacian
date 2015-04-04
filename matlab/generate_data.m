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

A=A+diag(-2*ones(1,m-1),-1)+diag(-2*ones(1,m-1),1)...
    +diag(-1*ones(1,m-2),2)+diag(-1*ones(1,m-2),-2);
%    diag(10*ones(1,m),0)+diag(-2*ones(1,m-1),1)...
A(1,1)=3;
A(2,2)=5;
A(3,3)=6;
A(4,4)=5;
A(5,5)=3;

for i=1:m
    for j=1:m
        fprintf(myFile,' %f ', A(i,j));
    end
    fprintf(myFile,'\n');
end
rang=rank(A)
b=[1; 0; 0; 0; -1];
A\b-0.5238*ones(5,1)






