function X=lyap3(A,B,C)
%Sylvester
[nr,nc]=size(C);
A0=kron(A,eye(nc))+kron(eye(nr),B');
try
    C1=C';
    X0=-inv(A0)*C1(:);
    X=transpose(reshape(X0,nc,nr));
catch
    error('Singular!');
end