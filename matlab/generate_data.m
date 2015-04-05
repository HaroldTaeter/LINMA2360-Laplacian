%
% Laplacian: 
% generate SDD matrix
%
clear all; 
close all; 
m=1000;

%% Print matrix A
file='dataL.txt';
myFile = fopen(file,'w+');
fprintf(myFile,'Number of nodes %d \n', m);
fprintf(myFile,'Number of edges %d \n', (m-1)+(m-2));

A=zeros(m,m);

A=A+diag(-2*ones(1,m-1),-1)+diag(-2*ones(1,m-1),1)...
    +diag(-1*ones(1,m-2),2)+diag(-1*ones(1,m-2),-2);

for i=1:m
    sum=0;
    for j=[1:i-1 , i+1:m]
        sum = sum + A(i,j);
    end   
    A(i,i)=-sum; 
end

% for i=1:m
%     for j=1:m
%         fprintf(myFile,' %f ', A(i,j));
%     end
%     fprintf(myFile,'\n');
% end

%% print vector b
% file2='datab.txt';
% myFile2 = fopen(file2,'w+');
% fprintf(myFile2,'%f \n',1);
% for i=1:m-2
% fprintf(myFile2,'%f \n',0);
% end
% fprintf(myFile2,'%f \n',-1);

%% Solve system in Matlab
rang=rank(A)
b=[1; zeros(m-2,1); -1];
tic();
v=A\b;
solution=v-v(1)*ones(m,1)
time_elapsed=toc()





