clear all; 
close all; 

m=9000; 
Ac=sparse(m,m);
Ac=sparse(Ac+diag(-2*ones(1,m-1),-1)+diag(-2*ones(1,m-1),1)...
    +diag(-1*ones(1,m-2),2)+diag(-1*ones(1,m-2),-2));

% PRINT
A=tril(Ac);
[row,col,val] = find(sparse(A));
fileID = fopen('edges9000.txt','w');
fprintf(fileID,'%d \n', m);
fprintf(fileID,'%d \n', (m-1)+(m-2));
fprintf(fileID, '%d %d %f\n', [row-ones(size(row)),col-ones(size(col)),-val]');
fclose(fileID);

% for i=1:m
%     Ac(i,i)=-sum(Ac(i,:));
% end

%% Print vector b
b=zeros(m,1);
b(1,1)=1;
b(m,1)=-1;
%file='datab1000.txt';
myFile = fopen('datab9000.txt','w+');
%for i=1:m
    fprintf(myFile,' %f \n', b);    
%end
fclose(myFile);

%% Solve system in Matlab
% tic();
% b=zeros(m,1);
% b(1,1)=1;
% b(m,1)=-1;
% Ac=sparse(Ac);
% sol=Ac\b;
% sol=sol-sol(1)*ones(m,1);
% time_elapsed=toc()


