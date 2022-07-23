function dy = classkcontrol_ljc_stable(t,x0,Inp,filename)
persistent gammashist
persistent time_counter
if nargin==3
% time_counter = 0;
time_counter = [time_counter 1];
Cs = Inp.Cs;Cb = Inp.Cb;W = Inp.W;
s0 = Inp.s0;b0 = Inp.b0;P=Inp.P;D=Inp.D;L=Inp.L;R=Inp.R;a=Inp.a;b=Inp.b;
tf = Inp.tf;
t_total = length(Inp.Yt);

t_span = tf/t_total;
gammashist(:,1)=zeros(size(Cs,1),1);
Yt = Inp.Yt;
if iscell(Yt)

    Yt = trans_cell_to_mat(Yt);
    [a_test,~,~] = catch_error(Yt(:,:,1),Cb,Cs);
    if sum(size(a_test) - size(Yt(:,:,1))) ~= 0
        [a_test,b] = size(Yt(:,:,1));
        Yt_new = zeros(b,a_test,length(Yt(1,1,:)));
        for i =1:length(Yt(1,1,:))
            [a_test,~,~] = catch_error(Yt(:,:,1),Cb,Cs);
            Yt_new(:,:,i) =a_test;
        end
        Yt = Yt_new;
%     else
%         Yt_new = Yt;
    end
end
t
Yt_temp=Yt;
for i =1:(t_total/t_span)
    if (t<t_span*i)||(t==t_span*i)
%         666
%         disp(i)
        Yt = Yt_temp(:,:,i);
        break
    else
        continue
    end
end
X2 = x0(1:numel(x0)/2);
X2d = x0(numel(x0)/2+1:end);
[U,Ez,V] = svd(P);
rankEz = rank(Ez);
U1 = U(:,1:rankEz);
U2 = U(:,rankEz+1:end);
E0 = Ez(1:rankEz,1:rankEz);

X2 = reshape(X2,3,size(U2,2));
X2d = reshape(X2d,3,size(U2,2));


X1 = D*V*E0^-1;

N = [X1 X2]*U';
Nd = [0*X1 X2d]*U';
l_hat= b0;
n = size(N,2); % number of nodes
nb = size(Cb,1); % number of bars
ns = size(Cs,1); % number of strings
m=1; % mass of each bar = 1 kg

% B=N*Cb';                                     %%%%%%%%%%%
% diag(B'*B)                                   %%%%%%%%%%%
%%  bar length correction
[N,Nd]=bar_length_correction_ClassK(N,Nd,Cb,P,D,l_hat);
X2 = N*U2;
X2d = Nd*U2;
%% calculate lagrange multiplier
B=N*Cb';Bd=Nd*Cb';S=N*Cs';Sd=Nd*Cs';
%       diag(B'*B) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% l_hat = sqrt(diag(diag(B'*B)));
l_hat = b0;
Cr=(1/2)*abs(Cb);
m_hat = m.*eye(nb);
M = (1/12)*Cb'*m_hat*Cb + Cr'*m_hat*Cr;
% Minv = 3*Cb'*m_hat^-1*Cb + 4*Cr'*m_hat^-1*Cr;
Minv = M^-1;
gamma_hat = diag(gammashist(:,end));
% kC = P'*Cb';kD = Cb*Minv*P;kE=P'*Minv*P;E = eye(3);
% kA = -S*gamma_hat*Cs*Minv*P+B*diag(diag(0.5*l_hat^(-2)*B'*(S*gamma_hat*Cs-W)*Cb'-(1/12)*l_hat^(-2)*(m*eye(nb))*(Bd'*Bd)))*Cb*Minv*P+W*Minv*P;
% Omg=0;
% for i=1:size(kC,2)
%     Omg=Omg+(1/(2*l_hat(i,i)^2))*kron(kC(:,i)',kron(B(:,i),(B(:,i)*kD(i,:))'));
% end
% Omg = (Omg - [kron(kE,E(1,:));kron(kE,E(2,:));kron(kE,E(3,:))])\[kA(1,:)';kA(2,:)';kA(3,:)']; 
% Omg = reshape(Omg,3,[])
Omg = calclagmat(N,Cb,B,Bd,Cs,gamma_hat,Minv,P,W,l_hat,m_hat);

%% control to find gamma
EYE = eye(size(Cb,1));
biglambda = []; tau = [];
for i=1:size(Cb,1)
    biglambda = [biglambda;(-1/(2*l_hat(i,i)^2))*B(:,i)'*S*diag(Cs*Cb'*EYE(:,i))];
    tau = [tau;(1/(2*l_hat(i,i)^2))*B(:,i)'*(W+Omg*P')*Cb'*EYE(:,i)+ (1/(12*l_hat(i,i)^2))*m*norm(Bd(:,i))^2];
end
Acon = a*eye(size(R,2));
Bcon = b*eye(size(R,2));
BTu = L*(W+Omg*P')*Minv*R + L*Nd*R*Acon + (L*N*R-Yt)*Bcon;
BIGTAU=[];meu=[];
EYE = eye(size(R,2));
for i=1:size(R,2)
BIGTAU = [BIGTAU;L*(S*diag(Cs*Minv*R*EYE(:,i))+B*diag(Cb*Minv*R*EYE(:,i))*biglambda)];
meu = [meu;BTu*EYE(:,i) - L*B*diag(Cb*Minv*R*EYE(:,i))*tau];
end
options = optimoptions('lsqlin','Display','off');
% 
[gammas,~,residual,exitflag] = lsqlin(BIGTAU,meu,[],[],[],[],zeros(ns,1),0.5*ones(ns,1),[],options);
% % % gammas = linprog(ones(ns,1),[],[],BIGTAU,meu,zeros(ns,1),[])
% % if size(gammas,1)~=size(Cs,1)
% %     gammas = lsqlin(BIGTAU,meu,[],[],[],[],zeros(ns,1),1*ones(ns,1),[],options)
% % end
gammashist = [gammashist gammas];


%% convert to reduce order dynamics and integrate
gamma_hat = diag(gammas);
lambda_c = diag(-biglambda*gammas-tau);
% lambda_c = diag(diag(0.5*l_hat^-2*B'*(S*gamma_hat*Cs-W-Omg*P')*Cb'-(1/12)*l_hat^-2*(m*eye(nb))*Bd'*Bd));
K = Cs'*gamma_hat*Cs - Cb'*lambda_c*Cb;
% K = pinv(N)*(W*M^-1*U1 + Omg*P'*M^-1*U1)*pinv(U1);
M_til = U2'*M*U2; 
K_til = U2'*K*U2;
W_til = W*U2-X1*U1'*K*U2;
X2dd = (W_til-X2*K_til)*M_til^-1;
dy = [reshape(X2d,numel(X2d),1);reshape(X2dd,numel(X2dd),1)];
else
    save gammashist.mat
end
end