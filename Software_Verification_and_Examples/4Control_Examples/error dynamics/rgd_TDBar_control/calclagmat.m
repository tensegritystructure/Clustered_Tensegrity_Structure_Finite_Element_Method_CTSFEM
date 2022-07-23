function Omg = calclagmat(N,Cb,B,Bd,Cs,gamma_hat,Minv,P,W,l_hat,m_hat)
C = Cb*P;
A=[];
for i=1:size(P,2)
    for j=1:3
        lambdaochat=[];
        for k=1:size(Cb,1)
            lambdaochat = [lambdaochat,-0.5*B(j,k)*C(k,i)*l_hat(k,k)^-2];
        end
        lambdaochat = diag(lambdaochat);

        phiij = N*Cb'*lambdaochat*Cb*Minv*P + BIGTHETA(j,i,P)*P'*Minv*P;
        A = [A,reshape(phiij,[],1)];
    end
end
lambdaoohat = diag(diag(0.5*l_hat^(-2)*B'*(N*Cs'*gamma_hat*Cs-W)*Cb'-(1/12)*l_hat^(-2)*m_hat*Bd'*Bd));
phi0 = (N*Cs'*gamma_hat*Cs-N*Cb'*lambdaoohat*Cb-W)*Minv*P;
B = reshape(phi0,[],1);
Omg = pinv(A) *B;
Omg = reshape(Omg,3,[]);
end