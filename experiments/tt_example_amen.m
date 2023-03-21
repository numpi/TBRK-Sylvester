%Comparation of TT Sylvester Solver  and amen_solve2 for
% diffusion problem with different values of d.

n=1024;

tol=1e-8;

fileID = fopen('tt_amen_symm.txt','w');
fprintf(fileID,'& d & residual & time (s) \n');
fclose(fileID);

for d=3:5
A=cell(1,d);
for i=1:d
A{i}=spdiags(ones(n,1) * [-1 2 -1], -1:1, n, n);
end
B=tt_random(n,d,2);

AA=tt_matrix(full(A{1}));
for i=2:d
    AA=tkron(AA,tt_eye(size(A{i},1)))+...
        tkron(tt_eye(AA.n),tt_matrix(full(A{i})));
end
AA=round(AA,1e-13);

options.m=zeros(1,d);
options.M=zeros(1,d);
for i=1:d
    options.m(i)=eigs(A{i},1,'smallestreal','Maxiterations',1e5);
    options.M(i)=eigs(A{i},1,'largestreal','Maxiterations',1e5);
end
options.maxit=30;
options.tol=tol;
options.real=true;
options.period=3;

options.poles="det2";
tic
 [X,resval] = TT_Sylvester_adaptive(A,B,options);
time=toc;
res=tt_residual(A,B,X);

fileID = fopen('tt_amen_symm.txt','a');
fprintf(fileID,'Krylov & $%i$ & $%.2e$ & $%.2f$ \n',d, res,time);
fclose(fileID);

tic
X=amen_solve2(AA,B,tol);
time=toc;
res=tt_residual(A,B,X);

fileID = fopen('tt_amen_symm.txt','a');
fprintf(fileID,'amen & $%i$ & $%.2e$ & $%.2f$ \n',d, res,time);
fclose(fileID);
end

%------------------------------------------------

function res=tt_residual(A,C,X)
Y=C;
for i=1:C.d
    Y=Y-tt_tensor_matrix_product(X,A{i},i);
end
Y=round(Y,1e-13);
res=norm(Y)/norm(C);
end


