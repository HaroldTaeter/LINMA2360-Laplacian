%
% Laplacian: 
% generate b vector
%

m=10;
b=zeros(m,1);
b(1)=1;
b(10)=-1;

file='dataB.txt';
myFile = fopen(file,'w+');

for i=1:m
    fprintf(myFile,' %d \n', b(i));    
end
