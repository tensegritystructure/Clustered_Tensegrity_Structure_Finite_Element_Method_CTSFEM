%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%roller-driven planer T-bar%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
%EXAMPLE
clc; clear all; close all;
% Global variable
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','Steel_string');
material{1}='linear_elastic'; % index for material properties: multielastic, plastic.
material{2}=0; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

% cross section design cofficient
thick=6e-3;        % thickness of hollow bar
hollow_solid=0;          % use hollow bar or solid bar in minimal mass design (1)hollow (0)solid
c_b=0.1;           % coefficient of safty of bars 0.5
c_s=0.1;           % coefficient of safty of strings 0.3

% static analysis set
substep=1;                                     %substep
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=1;             % save data or not (1) yes (0)no
savevideo=1;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
%% N C of the structure
% Manually specify node positions of double layer prism.
N=[1-0.5 0 0;1 1 0;1 1 0;1+0.5 0 0;1 -1 0]';    

% Manually specify connectivity indices.
 C_s_in = [1 3;3 4;4 5;5 1];  % This is indicating that string connection
 C_b_in = [1 4;2 5];  % Similarly, this is saying bar 1 connects node 1 to node 2,

 % Convert the above matrices into full connectivity matrices.
 C_b = tenseg_ind2C(C_b_in,N);%%
 C_s = tenseg_ind2C(C_s_in,N);
 C=[C_b;C_s];
 [ne,nn]=size(C);        % ne:No.of element;nn:No.of node

 % C_sta and C_end
 C_sta=C;
 C_sta(find(C==1))=0;
 C_sta=abs(C_sta);
 C_end=C;
 C_end(find(C==-1))=0;

 % Plot the structure to make sure it looks right
 tenseg_plot(N,C_b,C_s);
 title('Clustered D-bar in edges and nodes');

 %% unique node
 [N_unique,ia,ic]=unique(N','row','stable');
 I_u_temp=eye(numel(ic));
 I_unique0=zeros(numel(ic),numel(ia));
 for i=1:numel(ic)
     I_unique0(:,ic(i))=I_unique0(:,ic(i))+I_u_temp(:,i);
 end
 I_unique=kron(I_unique0,eye(3));
 %% Boundary constraints-transformation matrix of nodal coordinate
 pinned_X=[]; pinned_Y=[]; pinned_Z=(1:numel(ia))';
 [Ia_min,Ib_min,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,numel(ia));
 E_na=I_unique*Ia_min;
 E_nb=I_unique*Ib_min;
 %% sliding distance vector of pulley and transformation matrix
sld=zeros(nn,1);
% transformation matrix of sliding distance vector
free_sld=3; fixed_sld=setdiff(1:nn,free_sld);
E_s=eye(nn);
E_sa=E_s(:,free_sld);
E_sb=E_s(:,fixed_sld);
% transformation matrix of generalized coordinate
E_qa=blkdiag(E_na,E_sa);
%% Group/Clustered information 
%generate group index
gr={[3,4]};     % number of elements in one group
% gr=[];                     %if no group is used
Gp=tenseg_str_gp(gr,C);    %generate group matrix
S=Gp';                      % clustering matrix

fig_handle=figure
tenseg_plot_CTS(N,C,[1,2],S,fig_handle);
%% equilibrium matrix and SVD
%Calculate equilibrium matrix and member length
H=N*C';                     % element's direction matrix
l=sqrt(diag(H'*H));         % elements' length
Cell_H=mat2cell(H,3,ones(1,size(H,2)));          % transfer matrix H into a cell: Cell_H
A_2=kron(C',eye(3))*blkdiag(Cell_H{:})/diag(l);     % equilibrium matrix
A_RDT=E_qa'*[A_2;C_end'-C_sta'];

%SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_RDT);

%%  cross section information
% cross sectional area
A=1e-5*[100*[1;1];ones(4,1)];
% young's modulus
E=2.06e11*ones(ne,1);
% external force
f_ena=zeros(size(E_na,2),1);
f_esa=zeros(size(E_sa,2),1);

% member force design
index_gp=[1];
fd=-1e3;
E_nb=eye(ne);
e_d=E_nb(:,index_gp); % e_d is the matrix to select group of member with designed force
z=(e_d'*V2)\(e_d'*(fd-pinv(A_RDT)*[f_ena;f_esa]));   %self-stress coefficient
t=pinv(A_RDT)*[f_ena;f_esa]+V2*z;

% member rest length
l0=E.*A.*l./(t+E.*A);
%% tangent stiffness matrix
Kn=kron(C'*diag(l.\t)*C,eye(3));

K_T=[Kn+A_2*diag(E.*A./l0)*A_2'-A_2*diag(l.\t)*A_2',A_2*diag(E.*A./l0)*(C_end-C_sta);...
    (C_end-C_sta)'*diag(E.*A./l0)*A_2',(C_end-C_sta)'*diag(E.*A./l0)*(C_end-C_sta)];
K_Tg=blkdiag(Kn-A_2*diag(l.\t)*A_2',zeros(numel(sld)));%geometry stiffness
K_Te=[A_2*diag(E.*A./l0)*A_2',A_2*diag(E.*A./l0)*(C_end-C_sta);...
    (C_end-C_sta)'*diag(E.*A./l0)*A_2',(C_end-C_sta)'*diag(E.*A./l0)*(C_end-C_sta)];

K_T=0.5*(K_T+K_T');
[K_mode,D1] = eig(E_qa'*K_T*E_qa);         % eigenvalue of tangent stiffness matrix
k=diag(D1);   
[k_sort,E_nb]=sort(k);
K_mode_sort=K_mode(:,E_nb);
% plot mode
num_plt=4:9;
plot_mode_CTS2(K_mode_sort(1:8,:),k_sort,N,E_na,C,1:2,S,l,'tangent stiffness matrix',...
    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)','N/m',num_plt,0.2,saveimg,2);
% plot_mode_CTS2(K_mode_sort(1:8,:),k_sort,N,Ia,C,1:2,S,l,'tangent stiffness matrix',...
%     'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)','N/m',num_plt,0.3,saveimg,3);

K_mode_sort(:,1)'*E_qa'*K_T*E_qa*K_mode_sort(:,1)
K_mode_sort(:,1)'*E_qa'*K_Tg*E_qa*K_mode_sort(:,1)
K_mode_sort(:,1)'*E_qa'*K_Te*E_qa*K_mode_sort(:,1)

%% tangent stiffness of Clustered Tsg
t_c=S'\t;E_c=S'\E;A_c=S'\A;l0_c=S*l0;q=l.\t;A_2ac=E_na'*A_2*S';
[Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS3(E_na,C,S,t_c,E_na'*A_2,E_c,A_c,l0,l);

K_mode(:,1)'*Kt_aa*K_mode(:,1)      %eigenvalue
K_mode(:,1)'*Kg_aa*K_mode(:,1)      % eigenvalue of geometry stiffness
K_mode(:,1)'*Ke_aa*K_mode(:,1)      % eigenvalue of material stiffness

% [Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS2(E_na,C,q,A_2ac,E_c,A_c,l0_c);
% plot the mode shape of tangent stiffness matrix
num_plt=4:8;
% plot_mode_CTS2(K_mode,k,N,Ia,C,1:2,S,l,'tangent stiffness matrix',...
%     'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)',num_plt,0.2,saveimg,3);
plot_mode(K_mode,k,N,E_na,C_b,C_s,l,'natrual vibration',...
    'Order of Vibration Mode','Frequency (Hz)',num_plt,0.2,saveimg);

%% compare tangent stiffness of RDT and CTS

%% statics analysis
substep=8;

ind_w=[];w=[];
ind_db=[]; dnb0=[];
%ind_dl0_c=[1,2,3,4]'; dl0_c=[-400,-300,200,100]';
ind_dl0_c=[]'; dl0_c=[]';
% ind_dl0_c=[1,2,3]'; dl0_c=[-40,-30,10]';
[w_t,dnb_t,l0_ct,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,l0_c,b,gravity,[0;9.8;0],C,mass);
data.w_t=w_t;  % external force
data.dnb_t=dnb_t;% forced movement of pinned nodes
data.l0_t=l0_ct;% forced movement of pinned nodes
data.N=N_out{end};
data.substep=substep;    % substep

data_out3=static_solver_RDT(data);
t_t3=data_out3.t_out;          %member force in every step
n_t3=data_out3.n_out;          %nodal coordinate in every step
N_out3=data_out3.N_out;
t_c_t3=pinv(S')*t_t3;