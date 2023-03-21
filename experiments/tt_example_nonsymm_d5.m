
n=1024;
d=5;

h=1/(n-1);
vi=0.1;

w1=@(x)(1 + (x + 1).^2/4);
w2=@(y) (1+y)/2;

t=linspace(0,1,n);

Phi{1}=diag(w1(t));
Phi{2}=diag(w2(t));


T=cell(1,d);
for i=1:d
T{i}=linspace(0,1,n);
end

T1=spdiags(ones(n,1) * [-1 2 -1], -1:1, n, n);
T1=(1/h^2)*T1;

A=cell(1,d);
for i=1:d
A{i}=-vi*T1;
end

B1=diag(ones(1,n-1),1)-diag(ones(1,n-1),-1);
B1=(1/(2*h))*B1;

for i=1:2
    A{i}=vi*A{i}+Phi{i}*B1;
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
str="TT_nonsymm_d5_";
 [X,resval] = TT_Sylvester_adaptive(A,B,options);
 m=zeros(size(resval,1),1);
 for i=1:size(resval,1)
     m(i)=mean(resval(i,1:3));
 end
%semilogy(m, resval(:,end)', 'r-')
dlmwrite(str+"det.dat",[m,resval(:,end)], '\t');

options.poles="det2";
 [X,resval] = TT_Sylvester_adaptive(A,B,options);
%hold on
 m=zeros(size(resval,1),1);
 for i=1:size(resval,1)
     m(i)=mean(resval(i,1:3));
 end
%semilogy(m, resval(:,end)', 'b-')
dlmwrite(str+"det2.dat",[m,resval(:,end)], '\t');

options.poles="poly";
 [X,resval] = TT_Sylvester_adaptive(A,B,options);
% hold on
 m=zeros(size(resval,1),1);
  for i=1:size(resval,1)
     m(i)=mean(resval(i,1:3));
 end
%semilogy(m, resval(:,end)', 'k-')
dlmwrite(str+"poly.dat",[m,resval(:,end)], '\t');

options.poles="ext";
 [X,resval] = TT_Sylvester_adaptive(A,B,options);
% hold on
 m=zeros(size(resval,1),1);
  for i=1:size(resval,1)
     m(i)=mean(resval(i,1:3));
 end
%semilogy(m, resval(:,end)', 'g-')
dlmwrite(str+"ext.dat",[m,resval(:,end)], '\t');

%legend('det','det2','poly', 'ext')


