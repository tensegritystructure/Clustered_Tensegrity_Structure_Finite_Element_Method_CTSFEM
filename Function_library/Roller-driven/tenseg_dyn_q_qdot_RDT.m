function [Yd,t,l,q]=tenseg_dyn_q_qdot_RDT(time,Y,data_in)
% Input:
%   t: current time
%   Y0: current X, Xd values:   Y0=[X;Xd];
% Output:
%   Yd=[Xd,Xdd]

% global l l_c q E t f q q_d l0 A stress strain l0_c
C=data_in.C;
E_qa=data_in.E_qa;
E_qb=data_in.E_qb;
% S=data_in.S;
index_b=data_in.index_b;
index_s=data_in.index_s;
rho=data_in.rho;
% w0=data.w;
% dXb=data.dXb;
consti_data=data_in.consti_data;
ksi=data_in.ksi;
material=data_in.material;

A=data_in.A;
l0_t=data_in.l0_t;
q0=data_in.q;
w_t=data_in.w_t;           % external force
dqb_t=data_in.dqb_t;       % forced node displacement
dqb_d_t=data_in.dqb_d_t;    %velocity of forced moved node
dqb_dd_t=data_in.dqb_dd_t;    %velocity of forced moved node

% M=data_in.M;            %this should be renewed 
% D=data_in.D;            %this need renew
dt=data_in.dt;

nn=data_in.nn;
nq=numel(Y)/2;
qa=Y(1:nq,:);               %free node cooridnate
qa_d=Y(nq+1:end,:);         %free node velocity
%%
% tspan=0:dt:tf;
ind=floor(time/dt)+1;
% dnb=dnb_t(:,ind);
% nb_d=dnb_d_t(:,ind); %this is the velocity of fixed node
% nb_dd=dnb_dd_t(:,ind); %this is the acceleration of fixed node

% Get current pinned nodes displacement 
if ischar(dqb_t)
    run(dqb_t) % if dnb is a string, run that script
elseif size(dqb_t,2)==1
    dqb = dqb_t; % dnb can be constant
else
    dqb = dqb_t(:,ind); % or dnb can be time-varying
end

% Get current pinned nodes velocity 
if ischar(dqb_d_t)
    run(dqb_d_t) % if dnb is a string, run that script
elseif size(dqb_d_t,2)==1
    qb_d = dqb_d_t; % dnb can be constant
else
    qb_d = dqb_d_t(:,ind); % or dnb can be time-varying
end

% Get current pinned nodes acceleration 
if ischar(dqb_dd_t)
    run(dqb_dd_t) % if dnb is a string, run that script
elseif size(dqb_dd_t,2)==1
    qb_dd = dqb_dd_t; % dnb can be constant
else
    qb_dd = dqb_dd_t(:,ind); % or dnb can be time-varying
end

% Get current external forces
if ischar(w_t)
    run(w_t) % if w_t is a string, run that script
elseif size(w_t,2)==1
    w = w_t; % w_t can be constant
else
    w = w_t(:,ind); % or W can be time-varying
end

% Get current rest length
if ischar(l0_t)
    run(l0_t) % if w_t is a string, run that script
elseif size(l0_t,2)==1
    l0 = l0_t; % w_t can be constant
else
    l0 = l0_t(:,ind); % or W can be time-varying
end
%%   calculate stiffness matrix
% w=w_t(:,ind);
qb=E_qb\q0+dqb;
q=E_qa*qa+E_qb*qb;
q_d=E_qa*qa_d+E_qb*qb_d;
sld=q(3*nn+1:end);              % sliding vector
l0s=l0-C*sld;

N=reshape(q(1:3*nn),3,[]);
H=N*C';                     % element's direction matrix
Cell_H=mat2cell(H,3,ones(1,size(H,2)));  % transfer matrix H into a cell: Cell_H

n=q(1:3*nn);
l=sqrt(sum((N*C').^2))'; %bar length

strain=(l-l0s)./l0s;        %strain of member
[E,stress]=stress_strain(consti_data,index_b,index_s,strain,material);
E_sct=stress./strain;           %secant modulus
t=stress.*A;         %member force


Kn=kron(C'*diag(t./l)*C,eye(3));                      %stiffness matrix

%% calculate mass matrix

mass=rho.*A.*l0s;        % mass vector
M=tenseg_mass_matrix_RDT(mass,C,N,l); % generate mass matrix


% E_tts=S'*E_sct;     %Young's modulus TTS
% A_tts=S'*A;     % Cross sectional area TTS
% l0_tts=(f+E_tts.*A_tts).\E_tts.*A_tts.*l;   %rest length, TTS
% 
% l0_c=S*l0_tts;
% % l0_c=l0;    %this is specially used for plastic material
% mass=rho.*A_tts.*l0_tts;
% M=tenseg_mass_matrix(mass,C,0); % generate mass matrix
%% calculate damping matrix
A_2=kron(C',eye(3))*blkdiag(Cell_H{:})/diag(l);     % equilibrium matrix
d_c=2/sqrt(3)*sqrt(rho)*A.*E.^0.5;                  % cricital damping 

D=[A_2;C']*diag(ksi.*d_c)*[A_2;C']';
%% calculate accerlation
qa_dd=(E_qa'*M*E_qa)\(E_qa\(w-M*E_qb*qb_dd-[Kn*n;C'*t]-D*q_d));      %dynamic equation
Yd=[qa_d;qa_dd];
