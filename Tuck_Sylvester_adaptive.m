function [X,C,resval] = Tuck_Sylvester_adaptive(A,G,B,options)
%[X,resval] = Tuck_Sylvester_adaptive(A,G,options)
%
%  Rational Krylov method for multiterm Sylvester equations that ensures
%  the last pole is equal to infinity at each step, to check the residual.
%
%  The function solves the multiterm Sylvester equation
%
%     X x_1 A{1} + ... + X x_d A{d}  = G x_1 B{1} ... x_d B{d}
%
%  where the RHS is represented in Tucker format and it is assumed that 0 is not
%  contained in the sum of the field of values of the matrices A{i}.
%
%Input parameters:
%  A is a cell(1,d), where A{i} is square matrix of size n(i);
%
%  G is a d- dimensional tensor of size bs(1) x ... x bs(d)
%
%  B is a cell(1,d), where B{i} is UNITARY matrix of size n(i) x bs(i)
%
%  options.maxit = maximum number of iterations of tensor block rational
%  Krylov (default 100),
%
%  options.tol = max final accuracy (in terms of relative residual),
%  default 1e-10.
%
%  Stopping criterion:
%   ||X x_1 A{1} + ... + X x_d A{d} -(G x_1 B{1} ... x_d B{d})||_F
%  ---------------------------------------- -----------------------< tol
%              || G x_1 B{1} ... x_d B{d}||_F
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
%  X is a d-dimensional tensor;
%  S is a cell(1,d);
% the approximate solution is given by X x_1 C{1} ... x_d C{d}
%  resval   history of the vector
%   [number of iteration of TBRK, relative residual norm],
% computed each options.period iterations.
%
%
% References:
% [1] Casulli, A., Robol, L., An effcient block rational Krylov solver for 
%     Sylvester equations with adaptive pole selection
% [2] Casulli, A., Tensorized block rational Krylov methods for multiterm
%     Sylvester equations

if nargout>=3
    resval=[];
end

if nargin==3
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

c=options.period;

d=length(A);
V=cell(1,d);
K=cell(1,d);
H=cell(1,d);
Ap=cell(1,d);
resA=zeros(1,d);
k=zeros(1,d);

bs=zeros(1,d);
n=zeros(1,d);
for i=1:d
    [n(i),bs(i)]=size(B{i});
end

normB=norm(G(:));

for i=1:d
    k(i)=1;
    [V{i},K{i},H{i}] = rat_krylov(A{i}, B{i}, inf);
    Ap{i}=H{i}(1:bs(i),1:bs(i))/K{i}(1:bs(i),1:bs(i));
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
while mean(k)<=maxit
    k=k+ones(1,d);
    if options.real
        s=zeros(1,d);
        for i=1:d
            [V{i},K{i},H{i},Ap{i}] = Swapped_update_real(A{i},V{i},...
                K{i},H{i},np{i}, Ap{i});
            s(i)=length(np{i});
        end
        w=num2cell(size(G)+bs.*s);
        G(w{:})=0;
        k=k+s-1;
    else
        for i=1:d
            [V{i},K{i},H{i},Ap{i}] = Swapped_update(A{i},V{i},...
                K{i},H{i},np{i}, Ap{i});
        end
       w=num2cell(size(G)+bs.*s);
        G(w{:})=0;
    end
    if (c==options.period || max(k)>= maxit)
        if length(n)==2
            Y=sylvester(Ap{1},Ap{2}',G);
        else
            Y = laplace_merge(Ap, G );
        end
    end
    c
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
        ss=size(Y);
        for i=1:d
            order = [i, 1:i-1, i+1:d];
            Z = permute(Y,order);
            Z = H{i}(end-bs(i)+1:end, :) / K{i}(1:end-bs(i), :)*...
                reshape( Z, ss(i), prod( ss( order(2:d) ) ) );
            resA(i) = norm(Z,"fro");
        end

        %Frobenius norm
        resnrm = norm(resA,2)/normB;
        string=' %d';
        string=join(repmat(string, 1,d));
        string= join(['Iteration' , string , ', res =%e\n']);
        fprintf(string,k, resnrm);
        if resnrm<options.tol
            X=Y;
            C=cell(1,d);
            for i=1:d
                C{i}=V{i}(:,1:end-bs(i));
            end
        if nargout>=3
            resval=[resval;[k,resnrm]];
        end
            return
        end
        if nargout>=3
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
C=cell(1,d);
for i=1:d
    C{i}=V{i}(:,1:end-bs(i));
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

for i=length(s):-1:1
    [G,~]...
        =qr(K(end-(i-1)*bs-2*bs+1:end-(i-1)*bs,end-(i-1)*bs-bs+1:end-(i-1)*bs));
    G=G';
    V(:,end-(i-1)*bs-2*bs+1:end-(i-1)*bs)=V(:,end-(i-1)*bs-2*bs+1:end-(i-1)*bs)*G';
    K(end-(i-1)*bs-2*bs+1:end-(i-1)*bs,end-(i-1)*bs-2*bs+1:end)...
        =G*K(end-(i-1)*bs-2*bs+1:end-(i-1)*bs,end-(i-1)*bs-2*bs+1:end);
    H(end-(i-1)*bs-2*bs+1:end-(i-1)*bs,end-(i-1)*bs-2*bs+1:end)...
        =G*H(end-(i-1)*bs-2*bs+1:end-(i-1)*bs,end-(i-1)*bs-2*bs+1:end);

    [G,~] = qr(H(end-(i-1)*bs-bs+1:end-(i-1)*bs, end-(i-1)*bs-2*bs+1:end-(i-1)*bs).');
    G=G(:,end:-1:1);
    H(:, end-(i-1)*bs-2*bs+1:end-(i-1)*bs)=H(:, end-(i-1)*bs-2*bs+1:end-(i-1)*bs)*G;
    K(:, end-(i-1)*bs-2*bs+1:end-(i-1)*bs)=K(:, end-(i-1)*bs-2*bs+1:end-(i-1)*bs)*G;

    e=zeros(k-(i-1)*bs,bs);
    e(end-bs+1:end,:)=eye(bs);
    Ap(1:end+bs,end+1:end+bs)=...
        H(1:end-(i-1)*bs-bs,1:end-(i-1)*bs)*(K(1:end-(i-1)*bs-bs,1:end-(i-1)*bs)\e);
    Ap(end-bs+1:end,1:end)=...
        H(end-(i-1)*bs-2*bs+1:end-(i-1)*bs-bs,1:end-(i-1)*bs)/K(1:end-(i-1)*bs-bs,1:end-(i-1)*bs);
end
end

%----------------------------------------------------------------------------

function [X, eigvals] = laplace_merge( A, B, nmin )
% LAPLACE_MERGE( A, B ) computes the solution of the Laplace equation
%   X x_1 A{1} + ... + X x_d A{d} = B.
% LAPLACE_MERGE( A, B, nmin ) allows to specify the minimal block
% size.
d = size(A,2);
n=size(B);
Q=cell(1,d);
if nargin < 3
    if d == 3
        nmin = 26;
    elseif d == 4
        nmin = 18;
    else
        nmin = 14;
    end
end

eigvals=cell(1,d);

for i = 1:d
    [Q{i},A{i}] = schur(A{i},'complex');
    order = [i, 1:i-1, i+1:d];
    B = permute(B,order);
    B = Q{i}'*reshape( B, n(i), prod( n( order(2:d) ) ) );
    B = reshape( B, n(order) ); B = ipermute(B,order);
    eigvals{i}=diag(A{i});  
end

X = laplace_tri( A, B, nmin );

for i = 1:d
    order = [i, 1:i-1, i+1:d];
    X = permute(X,order);
    X = Q{i}*reshape( X, n(i), prod( n( order(2:d) ) ) );
    X = reshape( X, n(order) ); X = ipermute(X,order);
end

end

function X = laplace_tri(A, B, nmin)
    n = size(B); d = length(n);
    if n(1)*n(2) <= nmin^2
        if d == 3
           AA = kroneckerize(A{1},A{2});
           BB = reshape( B, n(1)*n(2), n(3) );
           % Pre-2018 Matlab command:
           % XX = sylvslv('C', AA, A{3}.', BB); 
           XX = matlab.internal.math.sylvester_tri( AA, conj( A{3} ), BB, 'I', 'I', 'trans');
           X = B;
           X(:) = XX;
        else
           A{1} = kroneckerize(A{1},A{2});
           A(2:d-1) = A(3:d);
           A = A(1:d-1);
           B = reshape( B, [ n(1)*n(2), n(3:end) ] );
           X = laplace_tri( A, B, nmin );
           X = reshape( X, n );
        end
        return
    end
    
    [~,i] = max(n);
    
    k = floor( n(i) / 2 ); if A{i}(k+1,k), k = k + 1; end
    
    Arec = A; Arec{i} = Arec{i}(k+1:n(i),k+1:n(i));
    X = zeros(size(B));
    
    for j=1:d, ind2{j}=':'; end, ind2{i}=k+1:n(i); ind1=ind2; ind1{i}=1:k;
    X2 = laplace_tri( Arec, B(ind2{:}), nmin );
    X(ind2{:})=X2;
    
    order = [i, 1:i-1, i+1:d]; 
    X2 = permute(X2,order);
    B1 = A{i}(1:k,k+1:n(i)) * reshape( X2, n(i)-k, prod( n( order(2:d) ) ) );
    B1 = reshape( B1, [ k, n(order(2:d)) ] );
    B1 = ipermute(B1,order);
    B(ind1{:})=B(ind1{:})-B1;
    
    Arec = A; Arec{i} = Arec{i}(1:k,1:k);
    X(ind1{:})=laplace_tri( Arec, B(ind1{:}), nmin);
end

function K = kroneckerize(A, B)
% Kroneckerize(A, B) returns kron(eye,A) + kron(B,eye)

[ma,na] = size(A);
[mb,nb] = size(B);

K = zeros(ma*mb,na*nb);
jj = 1-nb:0;

EA = eye(na);

for i = 1:ma
        ik = i:na:(na*mb-na+i);
        K(ik,ik) = B;
end

for i = 1:mb
        ik = 1+(i-1)*ma:i*ma;
        K(ik,ik) = K(ik,ik) + A;
end

end

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


