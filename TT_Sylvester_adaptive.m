function [X,resval] = TT_Sylvester_adaptive(A,G,options)
%[X,resval] = TT_Sylvester_adaptive(A,G,options)
%
%  Rational Krylov method for multiterm Sylvester equations that ensures
%  the last pole is equal to infinity at each step, to check the residual.
%
%  The function solves the multiterm Sylvester equation
%
%     X x_1 A{1} + ... + X x_d A{d}  = G
%
%  assuming that 0 is not contained in the sum of the field of values of the
%  matrices A{i}.
%
%Input parameters:
%  A, is a cell(1,d), where A{i} are square matrices;
%
%  C is a d- dimensional tensor in Tensor Train format as described in
%  the TT-toolbox [3]
%
%  options.maxit = maximum number of iterations of tensor block rational
%  Krylov (default 100),
%
%  options.tol = max final accuracy (in terms of relative residual),
%  default 1e-10.
%
%  Stopping criterion:
%       ||X x_1 A{1} + ... + X x_d A{d} -C||_F
%       ---------------------------------------- < tol
%                 || C ||_F
%
%  options.poles determines how to adaptively chose poles. The possible
%  choices are options.poles='det', options.poles='det2', options.poles='poly'
%  default is 'det'.
%  The pole selection algorithm are described in detail in [1,2].
%
%  options.m and option.M are d dimensional array such that the ith entry
%  contains lower and  upper bounds, respectively, for the real part of the
%  eigenvalues of the Kronecker sum of the matrices A{j} for j not equal to i;
%
%  options.real='true' runs the algorithm for A{i} real for each i;
%
%  options.period specifies the number of iterations after which the
%  residual is periodically computed (default options.period=5).
%
%
% Output parameters:
%
%  X  is the solution tensor in TT format
%  resval   history of the vector
%   [number of iteration of TBRK, relative residual norm],
% computed each options.period iterations.
%
%
%
% References:
% [1] Casulli, A., Robol, L., ...
% [2] Casulli, A.,....
% [3] Oseledets, I. V. "TT-Toolbox software; see https://github.com/oseledets."

if nargout==2
    resval=[];
end

if nargin==2
    options=[];
end

if isfield(options, 'tol')==0
    options.tol=1e-10;
end

if isfield(options, 'real')==0
    options.real=false;
end

if isfield(options, 'maxit')==0
    options.maxit=100;
end
maxit=options.maxit;


if isfield(options, 'poles')==0
    options.poles="det";
end

if isfield(options, 'period')==0
    options.period=5;
end

if isfield(options, 'TT_SVD')==0
    options.TT_SVD=false;
end

c=options.period;

d=G.d;
V=cell(1,d);
K=cell(1,d);
H=cell(1,d);
bs=zeros(1,d);
Ap=cell(1,d);
resA=zeros(1,d);
k=zeros(1,d);
B=cell(1,d);


normB=norm(G);

B{1}=reshape(G.core(G.ps(1):G.ps(2)-1),G.n(1),G.r(2));
bs(1)=G.r(2);
for i=2:d
    B{i}=reshape(G.core(G.ps(i):G.ps(i+1)-1),G.r(i),G.n(i),G.r(i+1));
    B{i} = permute(B{i},[2,1,3]);
    B{i}=reshape(B{i},G.n(i),G.r(i)*(G.r(i+1)));
    bs(i)=G.r(i)*G.r(i+1);
end

C=G;
for i=1:d
    k(i)=1;
    [V{i},K{i},H{i}] = rat_krylov(A{i}, B{i}, inf);
    Ap{i}=H{i}(1:bs(i),1:bs(i))/K{i}(1:bs(i),1:bs(i));
    C=tt_tensor_matrix_product(C,V{i}(:,1:bs(i))',i);
end

np=cell(1,d);
if options.poles=="poly"
    for i=1:d
        np{i}=Inf;
    end
elseif options.poles=="ext"
    for i=1:d
        np{i}=0;
    end
else
    if isfield(options, 'm')==0
        options.m=zeros(1,d);
        for i=1:d
            options.m(i)=eigs(A{i},1,'smallestreal');
        end
        summA=options.m*ones(d,1);
        for i=1:d
            options.m(i)=summA-options.m(i);
        end
    end
    m=options.m;
    if isfield(options, 'M')==0
        options.M=zeros(1,d);
        for i=1:d
            options.M(i)=eigs(A{i},1,'largestreal');
        end
        sumMA=options.M*ones(d,1);
        for i=1:d
            options.M(i)=sumMA-options.M(i);
        end
    end
    M=options.M;
    poles =cell(1,d);
    for j=1:d
        if abs(m(j))<abs(M(j))
            np{j}=-real(m(j));
            poles{j}=[-real(M(j)),np{j}];
        else
            np{j}=-real(M(j));
            poles{j}=[-real(m(j)),np{j}];
        end
    end
end
Y=C;
while max(k)<=maxit
    k=k+ones(1,d);
    if options.real
        s=zeros(1,d);
        for i=1:d
            [V{i},K{i},H{i},Ap{i}] = Swapped_update_real(A{i},V{i},...
                K{i},H{i},np{i}, Ap{i});
            s(i)=length(np{i});
        end
        for i=1:d
            C=tt_tensor_matrix_product(C,eye(C.n(i)+bs(i)*s(i),C.n(i)),i);
        end
        k=k+s-1;
    else
        for i=1:d
            [V{i},K{i},H{i},Ap{i}] = Swapped_update(A{i},V{i},...
                K{i},H{i},np{i}, Ap{i});
        end
        for i=1:d
            C=tt_tensor_matrix_product(C,eye(C.n(i)+bs(i)*s(i),C.n(i)),i);
        end
    end
    if (c==options.period || max(k)>= maxit)
        if options.TT_SVD
            [Y,resY] = tt_SVD_Sylv(Ap,C,options);
        else
            options.x0=Y;
            for i=1:d
                options.x0=tt_tensor_matrix_product(options.x0,eye(C.n(i),options.x0.n(i)),i);
            end
            [Y,resY] = tt_amen_small(Ap,C,options);
        end
    end

    if options.poles=="det"
        np = poles_selection_det(m,M,Ap,poles,bs,options);
        for j=1:d
            poles{j}=[poles{j},np{j}];
        end
    elseif options.poles=="det2"
        np = poles_selection_det2(m,M,Ap,poles,bs,options);
        for j=1:d
            poles{j}=[poles{j},np{j}];
        end
    elseif options.poles=="poly"
        for i=1:d
            np{i} = inf;
        end
    elseif options.poles=="ext"
        for i=1:d
            if np{i} == inf
                np{i} = 0;
            else
                np{i} = inf;
            end
        end
    end
    if(c==options.period || max(k)>= maxit)
        for i=1:d
            Z=tt_tensor_matrix_product(Y,...
                H{i}(end-bs(i)+1:end, :) / K{i}(1:end-bs(i), :),i);
            resA(i) = norm(Z);
        end

        %Frobenius norm
        resnrm = norm(resA,2)/normB;
        string=' %d';
        string=join(repmat(string, 1,d));
        string= join(['Iteration' , string , ', res =%e\n']);
        fprintf(string,k, resnrm);
        if resnrm<options.tol
            X=Y;
            for i=1:d
                X=tt_tensor_matrix_product(X,V{i}(:,1:end-bs(i)),i);
            end
            if nargout==2
                resval=[resval;[k,resnrm]];
            end
            return
        end
        if resnrm<resY
            X=Y;
            for i=1:d
                X=tt_tensor_matrix_product(X,V{i}(:,1:end-bs(i)),i);
            end
            fprintf(['The algorithm needs more memory to guarantee the requested residual.' ...
                '\n The output sathisfies res=%.2e \n'], resY);
            return
        end
        if nargout==2
            resval=[resval;[k,resnrm]];
        end
        c=0;
    else
        c=c+1;
    end
end
if  max(k)>= maxit
    string=' %d';
    string=join(repmat(string, 1,d));
    string= join(['Stopped at iteration' , string , '\n']);
    fprintf(string, k);
end
X=Y;
for i=1:d
    X=tt_tensor_matrix_product(X,V{i}(:,1:end-bs(i)),i);
end
end

%--------------------------------------------------------------------------

function [V,K,H,Ap] = Swapped_update(A,V,K,H,s, Ap)
%Step of rational Krylov method to add the pole s
% and swap poles ensuring that that the last pole
% is equal to infinity.

bs = size(H, 1) - size(H, 2);
param.extend=bs;

[V,K,H] = rat_krylov(A,V,K,H,s,param);
k=size(K,2);
[G,K(end-2*bs+1:end,end-bs+1:end)]=qr(K(end-2*bs+1:end,end-bs+1:end));
G=G';
V(:,end-2*bs+1:end)=V(:,end-2*bs+1:end)*G';
H(end-2*bs+1:end,end-2*bs+1:end)=G*H(end-2*bs+1:end,end-2*bs+1:end);

[G,~] = qr(H(end-bs+1:end, end-2*bs+1:end).');
G=G(:,end:-1:1);
H(:, end-2*bs+1:end)=H(:, end-2*bs+1:end)*G;
K(:, end-2*bs+1:end)=K(:, end-2*bs+1:end)*G;

e=zeros(k,bs);
e(end-bs+1:end,end-bs+1:end)=eye(bs);
Ap(1:end+bs,end+1:end+bs)=H(1:end-bs,1:end)*(K(1:end-bs,1:end)\e);
Ap(end-bs+1:end,1:end)=H(end-2*bs+1:end-bs,1:end)/K(1:end-bs,1:end);


end


function [V,K,H,Ap] = Swapped_update_real(A,V,K,H,s, Ap)
%Step of rational Krylov method to add the pole s (and its conjugate if s
%is not real), and swap poles ensuring that that the last pole
% is equal to infinity.

bs = size(H, 1) - size(H, 2);
param.extend=bs;
param.real=1;
[V,K,H] = rat_krylov(A,V,K,H,s,param);
k=size(K,2);
l=length(s);

[Q,R]=qr(K(end-(l+1)*bs+1:end,end-l*bs+1:end));
K(end-(l+1)*bs+1:end,end-l*bs+1:end)=R;
V(:,end-(l+1)*bs+1:end)=V(:,end-(l+1)*bs+1:end)*Q;
H(end-(l+1)*bs+1:end,end-(l+1)*bs+1:end)=...
    Q'*H(end-(l+1)*bs+1:end,end-(l+1)*bs+1:end);
[Q,~]=qr(H(end:-1:end-(l+1)*bs+1,end-(l+1)*bs+1:end)');
Q=Q(:,end:-1:1);
H(:,end-(l+1)*bs+1:end)=H(:,end-(l+1)*bs+1:end)*Q;
K(1:end-bs,end-(l+1)*bs+1:end)=...
    K(1:end-bs,end-(l+1)*bs+1:end)*Q;

e=zeros(k,l*bs);
e(end-l*bs+1:end,:)=eye(l*bs);
Ap(1:end+l*bs,end+1:end+l*bs)=...
    H(1:end-bs,1:end)*(K(1:end-bs,1:end)\e);
Ap(end-l*bs+1:end,1:end)=...
    H(end-(l+1)*bs+1:end-bs,1:end)/K(1:end-bs,1:end);
end
%------------------------------------

function  np = poles_selection_det(m,M,A,poles,bs,options)

d=length(A);
if isfield(options, 'eigvals')==0
    eigvals=cell(1,d);
    for j=1:d
    eigvals{j}=(eig(A{j}));
    end
else
    eigvals=options.eigvals;
end

for j=1:d
    eigvals{j}=sort(eigvals{j});
end
sumvals=0;
k=zeros(1,d);
eigB=cell(1,d);
for i=1:d
    k(i)=floor(size(A{i},1)/min(nthroot(1e3,d),size(A{i},1)));
    order = [i, 1:i-1, i+1:d];
    sumvals=sumvals+ipermute(eigvals{i}(1:k(i):end),order);
end

for i=1:d
    order = [i, 1:i-1, i+1:d];
    eigB{i}=permute(sumvals,order)-eigvals{i}(1:k(i):end);
    eigB{i}=eigB{i}(1,:);
end

np=cell(1,d);
for i=1:d
    np{i} = newpole_adaptive_det(-m(i),-M(i),...
        -eigB{i},eigvals{i},poles{i}, bs(i));
    if options.real
        if ~isreal(np{i}(1))
            np{i}(2)=conj(np{i}(1));
        end
    end
end
end

function np = newpole_adaptive_det(a,b,eigenvaluesA,eigenvaluesB, poles,bs)
%Computes the newpole for rational Krylov method maximizing the determinant.

poles = poles(:);

eigenvaluesB = eigenvaluesB(:);

poles = kron(ones(bs, 1), poles);

if isreal(eigenvaluesA)
    %t=linspace(a,b,40*length(eigenvaluesB));
    %t=linspace(a,b,500);
    %ma=min(abs(a),abs(b)); Ma=max(abs(a),abs(b));
    %t = sign(a)*logspace(log10(ma), log10(Ma), 100);
    %[~, jx] = max(abs(ratfun(t, eigenvaluesB, poles)));
    %np=t(jx);
    
    eHpoints=sort([a;b;eigenvaluesA(:)]);
    maxvals=zeros(2,length(eHpoints)-1);
    for j=1:length(eHpoints)-1
        t=linspace(eHpoints(j),eHpoints(j+1),20);
        [maxvals(1,j),jx]= max(abs(ratfun(t, eigenvaluesB, poles)));
        maxvals(2,j)=t(jx);
    end
    [~,jx]=max(maxvals(1,:));
    np=maxvals(2,jx);
else
    x=[eigenvaluesA(:);a;b];
    k = convhull(real(x),imag(x));
    maxvals=zeros(2,length(k)-1);
    for i=1:length(k)-1
        t=linspace(x(k(i)),x(k(i+1)),20);
        [maxvals(1,i),jx] =max(abs(ratfun(t, eigenvaluesB, poles)));
        maxvals(2,i)=t(jx);
    end
    [~,jx]=max(maxvals(1,:));
    np=maxvals(2,jx);
end
end

function r=ratfun(x,eH,s)

r=zeros(size(x));

for j=1:length(x)
    %r(j)=abs(prod( (x(j)-s)./(x(j)-eH) ));

    %to add infinity pole
    r(j)=abs(prod( (x(j)-s(2:end))./(x(j)-eH(2:end)) )/((x(j)-eH(1))));
end

return

end

function  np = poles_selection_det2(m,M,A,poles,bs,options)

d=length(A);
if isfield(options, 'eigvals')==0
    eigvals=cell(1,d);
    for j=1:d
    eigvals{j}=(eig(A{j}));
    end
else
    eigvals=options.eigvals;
end

for j=1:d
    eigvals{j}=sort(eigvals{j});
end
sumvals=0;
k=zeros(1,d);
eigB=cell(1,d);
for i=1:d
    k(i)=floor(size(A{i},1)/min(nthroot(1e3,d),size(A{i},1)));
    order = [i, 1:i-1, i+1:d];
    sumvals=sumvals+ipermute(eigvals{i}(1:k(i):end),order);
end

for i=1:d
    order = [i, 1:i-1, i+1:d];
    eigB{i}=permute(sumvals,order)-eigvals{i}(1:k(i):end);
    eigB{i}=eigB{i}(1,:);
end

np=cell(1,d);
for i=1:d
    np{i} = newpole_adaptive_det2(-m(i),-M(i),...
        -eigB{i},eigvals{i},poles{i}, bs(i));
    if options.real
        if ~isreal(np{i}(1))
            np{i}(2)=conj(np{i}(1));
        end
    end
end
end
%-----------------------------------------------------------------------
function np = newpole_adaptive_det2(a,b,eigenvaluesA, eigenvaluesB, poles,bs)
%Computes the newpole for rational Krylov method maximizing the product
%of a selected set of eigenvalues.

poles = poles(:);

eigenvaluesB = eigenvaluesB(:);

if isreal(eigenvaluesA)
    %t=linspace(a,b,40*length(eigenvaluesB));
    %ma=min(abs(a),abs(b)); Ma=max(abs(a),abs(b));
    %t = sign(a)*logspace(log10(ma), log10(Ma), 100);
    %     vals=zeros(1,length(t));
    %     for j=1:length(t)
    %         %[~,I]=sort(t(j)-eigenvaluesB, 'descend');
    %         %sorteig=eigenvaluesB(I);
    %         vals(j)=abs(prod( (t(j)-poles(2:end))./(t(j)-sorteig(bs+1:bs:end)) )...
    %             /((t(j)-sorteig(1))));
    %     end
    %     [~, jx] = max(vals);
    %     np=t(jx);

    eHpoints=sort([a;b;eigenvaluesA(:)]);
    maxvals=zeros(2,length(eHpoints)-1);
    for i=1:length(eHpoints)-1
        t=linspace(eHpoints(i),eHpoints(i+1),20);
        vals=zeros(1,length(t));
        for j=1:length(t)
            [~,I]=sort(abs(t(j)-eigenvaluesB), 'ascend');
            sorteig=eigenvaluesB(I);
            vals(j)=abs(prod( (t(j)-poles(2:end))./(t(j)-sorteig(bs+1:bs:end)) )...
                /((t(j)-sorteig(1))));
        end
        [maxvals(1,i),jx]= max(vals);
        maxvals(2,i)=t(jx);
    end
    [~,jx]=max(maxvals(1,:));
    np=maxvals(2,jx);

else
    x=[eigenvaluesA(:);a;b];
    k = convhull(real(x),imag(x));
    maxvals=zeros(2,length(k)-1);
    for i=1:length(k)-1
        t=linspace(x(k(i)),x(k(i+1)),20);
        vals=zeros(1,length(t));
        for j=1:length(t)
            [~,I]=sort(abs(t(j)-eigenvaluesB), 'ascend');
            sorteig=eigenvaluesB(I);
            vals(j)=abs(prod( (t(j)-poles(2:end))./(t(j)-sorteig(bs+1:bs:end)) )...
                /((t(j)-sorteig(1))));
        end
        [maxvals(1,i),jx] = max(vals);
        maxvals(2,i)=t(jx);
    end
    [~,jx]=max(maxvals(1,:));
    np=maxvals(2,jx);
end

end

function [X] = tt_tensor_matrix_product(X,A,i)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Z=reshape(X.core(X.ps(i):X.ps(i+1)-1),X.r(i),X.n(i),X.r(i+1));
Z=A*reshape(permute(Z,[2,1,3]),X.n(i),X.r(i)*X.r(i+1));
Z=ipermute(reshape(Z,size(A,1), X.r(i),X.r(i+1)),[2,1,3]);
X.core=[X.core(1:X.ps(i)-1);Z(:);X.core(X.ps(i+1):end)];
X.ps(i+1:end)=X.ps(i+1:end)+(X.ps(i)-X.ps(i+1)+length(Z(:)));
X.n(i)=size(A,1);
end

%----------------------------------------------------------------------

function [X,res] = tt_amen_small(A,C,options)
%[X,res] = tt_amen_small(A,C,options)
%
% The function solves the multiterm Sylvester equation
%
%     X x_1 A{1} + ... + X x_d A{d}  = C.
%
% It is indicated for matrices A{i} with moderate size (depending on the
% available memory). The method solves the full equation if there is 
% enough memory, otherwise uses the amen iteration.
%
%Input parameters:
%  A, is a cell(1,d), where A{i} are square matrices;
%
%  C is a d- dimensional tensor in Tensor Train format as described in
%  the TT-toolbox [3]
%
%  options.memory is the number of bytes of the available memory (if not 
%  specified automatically considers the full memory of the machine)
%
%  options.tol = max final accuracy (in terms of relative residual),
%  default 1e-10.
%
% Output parameters:
%
%  X  is the solution tensor in TT format
%  res is the relative norm of the residual
%
%       ||X x_1 A{1} + ... + X x_d A{d} -C||_F
%       --------------------------------------
%                 || C ||_F
%
% References:
% [1] Casulli, A., Robol, L., ...
% [2] Casulli, A.,....
% [3] Oseledets, I. V. "TT-Toolbox software; see https://github.com/oseledets."

if nargin==2
    options=[];
end

if isfield(options, 'tol')==0
    options.tol=1e-10;
end
   AA=tt_matrix(full(A{1}));
for i=2:C.d
    AA=tkron(AA,tt_eye(size(A{i},1)))+...
        tkron(tt_eye(AA.n),tt_matrix(full(A{i})));
end
AA=round(AA,1e-13);
X=amen_solve2(AA,C,options.tol,'x0',options.x0);
res=tt_residual(A,C,X);
end

%----------------------------------------------
function res=tt_residual(A,C,X)
Y=C;
for i=1:C.d
    Y=Y-tt_tensor_matrix_product(X,A{i},i);
end
Y=round(Y,1e-13);
res=norm(Y)/norm(C);
end


