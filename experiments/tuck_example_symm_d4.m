%Effect of changing poles in Tucker Sylvester Solver for
% diffusion problem with d=4.

n=1024;
d=4;

str="Diffusion d="+int2str(d)+" n="+int2str(n)+".mat";
load(str);
B=X;

options.m=zeros(1,d);
options.M=zeros(1,d);
for i=1:d
    options.m(i)=eigs(A{i},1,'smallestreal','Maxiterations',1e5);
    options.M(i)=eigs(A{i},1,'largestreal','Maxiterations',1e5);
end
options.maxit=20;
options.tol=1e-16;
options.real=true;
options.poles="det";
options.period=3;
str="Tuck_symm_d4_";
 [X,C,resval] = Tuck_Sylvester_adaptive(A,G,B,options);
%semilogy(resval(:,1)', resval(:,end)', 'r-')
dlmwrite(str+"det.dat",[resval(:,1),resval(:,end)], '\t');

options.poles="det2";
 [X,C,resval] = Tuck_Sylvester_adaptive(A,G,B,options);
 %hold on
%semilogy(resval(:,1)', resval(:,end)', 'b-')
dlmwrite(str+"det2.dat",[resval(:,1),resval(:,end)], '\t');

options.poles="poly";
 [X,C,resval] = Tuck_Sylvester_adaptive(A,G,B,options);
% hold on
%semilogy(resval(:,1)', resval(:,end)', 'k-')
dlmwrite(str+"poly.dat",[resval(:,1),resval(:,end)], '\t');

options.poles="ext";
 [X,C,resval] = Tuck_Sylvester_adaptive(A,G,B,options);
% hold on
%semilogy(resval(:,1)', resval(:,end)', 'g-')
dlmwrite(str+"ext.dat",[resval(:,1),resval(:,end)], '\t');

%legend('det','det2','poly', 'ext')


