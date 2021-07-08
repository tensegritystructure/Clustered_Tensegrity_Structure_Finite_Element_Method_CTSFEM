function Yd=tenseg_control_ckclose_CTS(t,Y,data_in)
% Input:
%   t: current time
%   Y0: current X, Xd values:   Y0=[X;Xd];
% Output:
%   Yd=[Xd,Xdd]

global l q E f f_c n n_d l0 A stress strain M l0_c exitflag
C=data_in.C;
Ia=data_in.Ia;
Ib=data_in.Ib;
S=data_in.S;
index_b=data_in.index_b;
index_s=data_in.index_s;
rho=data_in.rho;
% w0=data.w;
% dXb=data.dXb;
consti_data=data_in.consti_data;
% slack=data_in.slack;
% plastic=data_in.plastic;
% multielastic=data_in.multielastic;
material=data_in.material;
lb=data_in.lb;      %lower bound
ub=data_in.ub;      %upper bound

A=data_in.A;
E=data_in.E;
l0=data_in.l0_t;           % rest length
n0=data_in.N(:);
w_t=data_in.w_t;           % external force
dnb_t=data_in.dnb_t;       % forced node displacement
dnb_d_t=data_in.dnb_d_t;    %velocity of forced moved node
dnb_dd_t=data_in.dnb_dd_t;    %velocity of forced moved node

D=data_in.D;
dt=data_in.dt;

% control parameters
Ic=data_in.Ic; n_ct_t=data_in.n_ct_t;n_ct_dt=data_in.n_ct_dt;n_ct_ddt=data_in.n_ct_ddt;
a_c=data_in.a_c;b_c=data_in.b_c;
I_act=data_in.I_act;I_pas=data_in.I_pas;

nf=numel(Y)/2;
na=Y(1:nf,:);               %free node cooridnate
na_d=Y(nf+1:end,:);         %free node velocity
%%
% tspan=0:dt:tf;
ind=floor(t/dt)+1;
% dnb=dnb_t(:,ind);
% nb_d=dnb_d_t(:,ind); %this is the velocity of fixed node
% nb_dd=dnb_dd_t(:,ind); %this is the acceleration of fixed node

% Get current pinned nodes displacement 
if ischar(dnb_t)
    run(dnb_t) % if dnb is a string, run that script
elseif size(dnb_t,2)==1
    dnb = dnb_t; % dnb can be constant
else
    dnb = dnb_t(:,ind); % or dnb can be time-varying
end

% Get current pinned nodes velocity 
if ischar(dnb_d_t)
    run(dnb_d_t) % if dnb is a string, run that script
elseif size(dnb_d_t,2)==1
    nb_d = dnb_d_t; % dnb can be constant
else
    nb_d = dnb_d_t(:,ind); % or dnb can be time-varying
end

% Get current pinned nodes acceleration 
if ischar(dnb_dd_t)
    run(dnb_dd_t) % if dnb is a string, run that script
elseif size(dnb_dd_t,2)==1
    nb_dd = dnb_dd_t; % dnb can be constant
else
    nb_dd = dnb_dd_t(:,ind); % or dnb can be time-varying
end

% Get current external forces
if ischar(w_t)
    run(w_t) % if w_t is a string, run that script
elseif size(w_t,2)==1
    w = w_t; % w_t can be constant
else
    w = w_t(:,ind); % or W can be time-varying
end

% % Get current rest length
% if ischar(l0_t)
%     run(l0_t) % if w_t is a string, run that script
% elseif size(l0_t,2)==1
%     l0 = l0_t; % w_t can be constant
% else
%     l0 = l0_t(:,ind); % or W can be time-varying
% end

% Get current control target
if  size(n_ct_ddt,2)==1
    n_ct = n_ct_t;%  constant   
    n_ct_d = n_ct_dt;
    n_ct_dd = n_ct_ddt;  
else
    n_ct = n_ct_t(:,ind); %  time-varying
    n_ct_d = n_ct_dt(:,ind);
    n_ct_dd = n_ct_ddt(:,ind);
end
%%   calculate stiffness matrix
% w=w_t(:,ind);
nb=Ib'*n0+dnb;
n=Ia*na+Ib*nb;
n_d=Ia*na_d+Ib*nb_d;

l=sqrt(sum((reshape(n,3,[])*C').^2))'; %bar length
l_c=S*l;

M_aa=Ia'*M*Ia;
H=reshape(n,3,[])*C';                     % element's direction matrix
Cell_H=mat2cell(H,3,ones(1,size(H,2)));          % transfer matrix H into a cell: Cell_H
A_2c=kron(C',eye(3))*blkdiag(Cell_H{:})*(diag(l)\S');     % equilibrium matrix

TAU=Ic'*(M_aa\Ia')*A_2c;
TAU_pas=TAU*I_pas;
TAU_act=TAU*I_act;
meu=Ic'*(M_aa\Ia')*(w-M*Ib*nb_dd-D*n_d)-n_ct_dd+a_c*(Ic'*na_d-n_ct_d)+b_c*(Ic'*na-n_ct);
t_pas=I_pas'*(E.*A./l0.*(l_c-l0));
% q_pas=(I_pas'*E).*(I_pas'*A)./(I_pas'*l0).*(I_pas'*l-I_pas'*l0);

 options = optimoptions('lsqlin','Display','off');
[t_act,~,residual,exitflag] = lsqlin(TAU_act,meu-TAU_pas*t_pas,[],[],[],[],lb,ub,[],options);%0*ones(size(I_act,2),1)
% max(residual)       % maximum residual
% t_act=pinv(TAU_act)*(meu-TAU_pas*t_pas);      % directly solve the t_act    

% if exitflag~=1
%     t_act=I_act'*f_c;   %use the value in last step
% end

f_c=I_pas*t_pas+I_act*t_act;       %the force vector
f=S'*f_c;
q_c=f_c./l_c;
q=f./l;

K=kron(C'*diag(q)*C,eye(3));                      %stiffness matrix

% strain=(l_c-l0)./l0;        %strain of members
% [E,stress]=stress_strain(consti_data,index_b,index_s,strain,material);
% f_c=stress.*A;         %member force
% f=S'*f_c;
% q_c=f_c./l_c;
% q=f./l;      %reculate force density
% q_bar=diag(q);
% K=kron(C'*q_bar*C,eye(3));                      %stiffness matrix

%% calculate mass matrix
E_tts=S'*E;     %Young's modulus TTS
A_tts=S'*A;     % Cross sectional area TTS
l0_tts=(f+E_tts.*A_tts).\(E_tts.*A_tts.*l);   %rest length, TTS
l0_c=S*l0_tts;
strain=(l_c-l0_c)./l0_c;
stress=f_c./A;
mass=rho.*A_tts.*l0_tts;
M=tenseg_mass_matrix(mass,C,0); % generate mass matrix

%% calculate accerlation
na_dd=(Ia'*M*Ia)\(Ia'*(w-M*Ib*nb_dd-K*n-D*n_d));      %dynamic equation
Yd=[na_d;na_dd];
