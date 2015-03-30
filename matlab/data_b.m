%
% Laplacian: 
% generate b vector
%

m=5;
b=zeros(m,1);
b(1)=1;
b(m)=-1;

file='datab.txt';
myFile = fopen(file,'w+');

for i=1:m
    fprintf(myFile,' %d \n', b(i));    
end
