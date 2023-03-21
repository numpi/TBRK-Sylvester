n=1024;
fileID = fopen('tt_increasing_d.txt','w');
fprintf(fileID,'$d$ & residual  &   time (s) \n');
fclose(fileID);
for d=5:5:20

A=cell(1,d);
for i=1:d
A{i}=spdiags(ones(n,1) * [-1 2 -1], -1:1, n, n);
end
B=tt_random(n,d,2);

options.m=zeros(1,d);
options.M=zeros(1,d);
for i=1:d
    options.m(i)=eigs(A{i},1,'smallestreal','Maxiterations',1e5);
    options.M(i)=eigs(A{i},1,'largestreal','Maxiterations',1e5);
end
options.maxit=50;
options.tol=1e-6;
options.real=true;
options.poles="det";
options.period=3;

tic
 [X,res] = TT_Sylvester_adaptive(A,B,options);
time=toc;

fileID = fopen('tt_increasing_d.txt','a');
fprintf(fileID,'$%i$ & $%.2e$ &  $%.2f$ \n',d,res(end),time);
fclose(fileID);
end




