%Time and number of iterations changing poles in Tucker Sylvester Solver for
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
options.maxit=100;
options.tol=1e-4;
options.real=true;
options.poles="det";
options.period=1;

tic
 [X,C,resval] = Tuck_Sylvester_adaptive(A,G,B,options);
 time_det_e4=toc;
 it_det_e4=resval(end,1:3);
 res_det_e4=resval(end,4);

options.poles="det2";
tic
 [X,C,resval] = Tuck_Sylvester_adaptive(A,G,B,options);
time_det2_e4=toc;
 it_det2_e4=resval(end,1:3);
 res_det2_e4=resval(end,4);

options.poles="ext";
tic
 [X,C,resval] = Tuck_Sylvester_adaptive(A,G,B,options);
time_ext_e4=toc;
 it_ext_e4=resval(end,1:3);
 res_ext_e4=resval(end,4);

 options.tol=1e-6;

 options.poles="det";
 tic
 [X,C,resval] = Tuck_Sylvester_adaptive(A,G,B,options);
 time_det_e6=toc;
 it_det_e6=resval(end,1:3);
 res_det_e6=resval(end,4);

options.poles="det2";
tic
 [X,C,resval] = Tuck_Sylvester_adaptive(A,G,B,options);
time_det2_e6=toc;
 it_det2_e6=resval(end,1:3);
 res_det2_e6=resval(end,4);

options.poles="ext";
tic
 [X,C,resval] = Tuck_Sylvester_adaptive(A,G,B,options);
time_ext_e6=toc;
 it_ext_e6=resval(end,1:3);
 res_ext_e6=resval(end,4);


disp("poles & iterations & time (s) & residual \\" )
 disp("det & $ " +int2str(it_det_e4)+" $ & $" + num2str(time_det_e4,'%.2f')+" $&$ " + ...
    num2str(res_det_e4,'%.2e')+" $ \\")
 disp("det2 & $" +int2str(it_det2_e4)+" $ &$ "+ num2str(time_det2_e4,'%.2f')+" $&$ " + ...
    num2str(res_det2_e4,'%.2e')+"$ \\")
  disp("ext & $ " +int2str(it_ext_e4)+" $ &$ "+ num2str(time_ext_e4,'%.2f')+" $&$ " + ...
    num2str(res_ext_e4,'%.2e')+" $ \\")
  disp(' ')
 disp("poles & iterations & time (s) & residual \\" )
  disp("det & $ " +int2str(it_det_e6)+" $ & $" + num2str(time_det_e6,'%.2f')+" $&$ " + ...
    num2str(res_det_e6,'%.2e')+" $ \\")
 disp("det2 & $" +int2str(it_det2_e6)+" $ &$ "+ num2str(time_det2_e6,'%.2f')+" $&$ " + ...
    num2str(res_det2_e6,'%.2e')+"$ \\")
  disp("ext & $ " +int2str(it_ext_e6)+" $ &$ "+ num2str(time_ext_e6,'%.2f')+" $&$ " + ...
    num2str(res_ext_e6,'%.2e')+" $ \\")






