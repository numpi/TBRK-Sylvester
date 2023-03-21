%Effect of changing poles in Tucker Sylvester Solver for
% convection-diffusion problem with d=3.

n=1024;
d=3;

str="Convection-Diffusion d="+int2str(d)+" n="+int2str(n)+".mat";
load(str);
B=X;
% 
% for i=1:d
% [B{i},~]=qr(randn(n,size(G,i)),0);
% end

%options.m=zeros(1,d);
%options.M=zeros(1,d);
m=zeros(d,1);
M=zeros(d,1);
for i=1:d
    e=eig(A{i});
    [~,I]=min(real(e));
    m(i)=e(I);
    [~,I]=max(real(e));
    M(i)=e(I);
    %options.m(i)=eigs(A{i},1,'smallestreal','Maxiterations',1e5);
    %options.M(i)=eigs(A{i},1,'largestreal','Maxiterations',1e5);
end
options.m=m;
options.M=M;
options.maxit=20;
options.tol=1e-16;
options.real=true;
options.poles="det";
options.period=3;
str="Tuck_nonsymm_d3_";
 [X,C,resval] = Tuck_Sylvester_adaptive(A,G,B,options);
 m=zeros(size(resval,1),1);
 for i=1:size(resval,1)
     m(i)=mean(resval(i,1:3));
 end
%semilogy(m, resval(:,end)', 'r-')
dlmwrite(str+"det.dat",[m,resval(:,end)], '\t');

options.poles="det2";
 [X,C,resval] = Tuck_Sylvester_adaptive(A,G,B,options);
 %hold on
 m=zeros(size(resval,1),1);
 for i=1:size(resval,1)
     m(i)=mean(resval(i,1:3));
 end
%semilogy(m, resval(:,end)', 'b-')
dlmwrite(str+"det2.dat",[m,resval(:,end)], '\t');

options.poles="poly";
 [X,C,resval] = Tuck_Sylvester_adaptive(A,G,B,options);
 %hold on
 m=zeros(size(resval,1),1);
  for i=1:size(resval,1)
     m(i)=mean(resval(i,1:3));
 end
%semilogy(m, resval(:,end)', 'k-')
dlmwrite(str+"poly.dat",[m,resval(:,end)], '\t');

options.poles="ext";
 [X,C,resval] = Tuck_Sylvester_adaptive(A,G,B,options);
% hold on
 m=zeros(size(resval,1),1);
  for i=1:size(resval,1)
     m(i)=mean(resval(i,1:3));
 end
%semilogy(m, resval(:,end)', 'g-')
dlmwrite(str+"ext.dat",[m,resval(:,end)], '\t');
%legend('det','det2','poly', 'ext')


