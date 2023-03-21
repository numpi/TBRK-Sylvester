%Effect of changing poles in TT Sylvester Solver for
% diffusion problem with d=6.

n=1024;
d=6;


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
options.maxit=30;
options.tol=1e-10;
options.real=true;
options.poles="det";
options.period=3;

str="TT_symm_d6_";
 [X,resval] = TT_Sylvester_adaptive(A,B,options);
%semilogy(resval(:,1)', resval(:,end)', 'r-')
dlmwrite(str+"det.dat",[resval(:,1),resval(:,end)], '\t');

options.poles="det2";
 [X,resval] = TT_Sylvester_adaptive(A,B,options);
%hold on
%semilogy(resval(:,1)', resval(:,end)', 'b-')
dlmwrite(str+"det2.dat",[resval(:,1),resval(:,end)], '\t');

options.poles="poly";
 [X,resval] = TT_Sylvester_adaptive(A,B,options);
% hold on
%semilogy(resval(:,1)', resval(:,end)', 'k-')
dlmwrite(str+"poly.dat",[resval(:,1),resval(:,end)], '\t');

options.poles="ext";
 [X,resval] = TT_Sylvester_adaptive(A,B,options);
% hold on
%semilogy(resval(:,1)', resval(:,end)', 'g-')
dlmwrite(str+"ext.dat",[resval(:,1),resval(:,end)], '\t');

%legend('det','det2','poly', 'ext')


