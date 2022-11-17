%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% origami+D bar tensegrity with bottom bracing

%EXAMPLE
clc; clear all; close all;
% Global variable
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','Steel_string');
material{1}='linear_elastic'; % index for material properties: multielastic, plastic.
material{2}=0; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

% cross section design cofficient
thick=6e-4;        % thickness of hollow bar
hollow_solid=0;          % use hollow bar or solid bar in minimal mass design (1)hollow (0)solid
c_b=0.1;           % coefficient of safty of bars 0.5
c_s=0.1;           % coefficient of safty of strings 0.3

% static analysis set
substep=15;                                     %荷载子步
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=1;             % save data or not (1) yes (0)no
savevideo=1;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
% move_ground=0;          % for earthquake, use pinned nodes motion(1) or add inertia force in free node(0) 
%% %% N C of the structure
% Manually specify node positions
% N=[0 1 0;1 1 0;0 0 0;1 0 0;0.7 0.3 0]';
p=3;        %complexity of shelter
b1=0.2;b2=0.2;h=0.1;        % geometry coefficient
N=zeros(3,8*p+6);           % initialize N
N(1,1:2:2*p+1)=[b1*(0:p)];
N(:,2:2:2*p)=[b1*(0.5+0:p);zeros(1,p);h*ones(1,p)];
N(:,2*p+2:4*p+2)=N(:,1:2*p+1)+[0;b2;0]*ones(1,2*p+1);
N(:,4*p+3:5*p+2)=diag([1,1,-1])*N(:,2:2:2*p);
N(:,5*p+3:6*p+2)=diag([1,1,-1])*N(:,2*p+3:2:4*p+1);
%center node
N(:,6*p+3:8*p+2)=[b1/2*(0.5:2*p-0.5);b2/2*ones(1,2*p);h/2*ones(1,2*p)];
%pinned node
N(:,8*p+3:8*p+6)=[0 0 norm([h,b1/2]);0 0 -norm([h,b1/2]);0 b2 norm([h,b1/2]);0 b2 -norm([h,b1/2])]';


% Manually specify connectivity indices.
% C_b_in_v=[1 2;2 4;4 3;3 1];% vertical horizontal bar
% C_b_in_d=[1 5;2 5;4 5;3 5];% diagonal bar
% C_in = [C_b_in_v;C_b_in_d];   % This is indicating the bar connection
C_b_in_1=[[1:2*p,2*p+2:4*p+1];[2:2*p+1,2*p+3:4*p+2]]';% bars in top
C_b_in_2=[[1:2:2*p-1,3:2:2*p+1];kron(ones(1,2),4*p+3:5*p+2)]'; % bottom bar
C_b_in_3=[[2*p+2:2:4*p,2*p+4:2:4*p+2];kron(ones(1,2),5*p+3:6*p+2)]'; % bottom bar 2
C_b_in_4=[1 2*p+1;2*p+2 4*p+2]';     %2 edge 
C_b_in_5=[2:2*p;2*p+3:4*p+1]';  % rotational hinge
% bending hinge
C_b_in_6=kron(ones(2*p,1),[1 2 2*p+2 2*p+3;(6*p+3)*ones(1,4)]')+kron([0:2*p-1]',ones(4,2));
 % diagonal braces in bottom
C_b_in_7=[1:2:2*p-1,4*p+3:5*p+2,2*p+2:2:4*p,5*p+3:6*p+2;5*p+3:6*p+2,2*p+4:2:4*p+2,4*p+3:5*p+2,3:2:2*p+1]';% diagonal braces in bottom
% strings
C_s_in_1=[8*p+3 2:2:2*p-2 8*p+5 2*p+3:2:4*p-1;2:2:2*p 2*p+3:2:4*p+1]'; %top string
C_s_in_2=[8*p+4 4*p+3:5*p+1 8*p+6 5*p+3:6*p+1;4*p+3:5*p+2 5*p+3:6*p+2]'; %bottom string
C_s_in_3=[[2:2:2*p 2*p+3:2:4*p+1];[4*p+3:6*p+2]]';                      %vertical string
% C_in = [C_b_in_1;C_b_in_2;C_b_in_3;C_b_in_4;C_b_in_5;C_b_in_6;C_s_in_1;C_s_in_2];
C_b_in=[C_b_in_1;C_b_in_2;C_b_in_3;C_b_in_4;C_b_in_5;C_b_in_6;C_b_in_7];
C_s_in=[C_s_in_1;C_s_in_2;C_s_in_3];
C_in=[C_b_in;C_s_in];
% Convert the above matrices into full connectivity matrices.
C = tenseg_ind2C(C_in,N);
[ne,nn]=size(C);        % ne:No.of element;nn:No.of node

C_b=tenseg_ind2C(C_b_in,N);C_s=tenseg_ind2C(C_s_in,N);
n_b=size(C_b,1);n_s=size(C_s,1);        % number of bars and string
% connectivity matrix for plot
C_bar_in=[C_b_in_1;C_b_in_2;C_b_in_3;C_b_in_7];      % real bars
C_bar=tenseg_ind2C(C_bar_in,N);
C_rot_h_in=[C_b_in_4;C_b_in_5];             %rotational hinges
C_rot_h=tenseg_ind2C(C_rot_h_in,N);

% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,C_s);


%% define hinge, rigid hinge
% C_in_h is the connectivity of higes, can be written in a function!!!!!!!!!
C_in_h=[C_b_in_5;C_b_in_6];
n_h=size(C_in_h,1);         % number of hinge

[~,index_h]=ismember(C_in_h,C_in,'rows');   % index_h is the index number of hinge
% index_rh=[];index_rh_in_h=[];
[~,index_rig_h]=ismember(C_b_in_6,C_in,'rows');   % index_h is the index number of rigid hinge
[~,index_rot_h]=ismember(C_b_in_5,C_in,'rows');   % index_h is the index number of rotational hinge

[~,index_roth_in_h]=ismember(C_b_in_5,C_in_h,'rows');   % index_h is the index of rigid hinge in all hinge
[~,index_righ_in_h]=ismember(C_b_in_6,C_in_h,'rows');   % index_h is the index of rigid hinge in all hinge

C_h=tenseg_ind2C(C_in(index_h,:),N);     % connectivity matrix of all edges
C_rig_h=tenseg_ind2C(C_in(index_rig_h,:),N); % connectivity matrix of rigid edges
% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,C_s);


%% connectivity of triangle element Ca
% Ca can be written in a function!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Ca=kron(ones(1,2*p),[6*p+3 1 2;6*p+3 2 2*p+3;6*p+3 2*p+3 2*p+2;;6*p+3 2*p+2 1]')+kron([0:2*p-1],ones(3,4));% Ca=generate_Ca(C_in,N);
% Ca=zeros(3,1)

[~,np]=size(Ca);        % ne:No.of element;np:No.of plate

% plot the origami configuration
tenseg_plot_ori(N,C_bar,C_s,C_rot_h,C_rig_h,[],[],[0,30],[] ,[],Ca);
axis off
% plot the origami configuration without bars and strings
tenseg_plot_ori(N,[],[],C_rot_h,C_rig_h,[],[],[0,30],[] ,[],Ca);
axis off
% plot the bars and strings
tenseg_plot_ori(N,C_bar,C_s,[],[],[],[],[0,30],[] ,[],[]);
ylim([-0.1,0.1])
axis off
%% transformation matrix from element to structure

E_n=cell(1,n_h);            %transformation matrix from element node to total node
node_in_hinge=zeros(n_h,4);
I=eye(3*nn);

for i=1:n_h
node2=C_in_h(i,1);  % start node of the hinge
node3=C_in_h(i,2);  % end node of the hinge
for j=1:np
    if (node2==Ca(1,j)&node3==Ca(2,j))|(node2==Ca(2,j)&node3==Ca(3,j))|(node2==Ca(3,j)&node3==Ca(1,j))
        node1=setdiff(Ca(:,j),[node2;node3]);
    elseif (node2==Ca(2,j)&node3==Ca(1,j))|(node2==Ca(3,j)&node3==Ca(2,j))|(node2==Ca(1,j)&node3==Ca(3,j))
        node4=setdiff(Ca(:,j),[node2;node3]);
    end
end
node_in_hinge(i,:)=[node1,node2,node3,node4];
E_n{i}=I(:,kron(node_in_hinge(i,:),3*ones(1,3))-kron(ones(1,4),[2,1,0]));
end
E_n_total=cell2mat(E_n);        % transformation matrix of the whole structure
%% Boundary constraints
pinned_X=[1 2*p+2 8*p+3:8*p+6]'; pinned_Y=[1 2*p+2 8*p+3:8*p+6]'; pinned_Z=[1 2*p+2 8*p+3:8*p+6]';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);

%% generate group index for tensegrity torus structure
gr=[];                      %no group is used
Gp=tenseg_str_gp(gr,C);    %generate group matrix
S=Gp';                      % clustering matrix
%% equilibrium matrix

% equilibrium matrix of truss
[A_1,A_1c,A_1a,A_1ac,A_2,A_2c,A_2a,A_2ac,l,l_c]=tenseg_equilibrium_matrix_CTS(N,C,S,Ia);
% [A_1,A_1g,A_2,A_2g,l,l_gp]=tenseg_equilibrium_matrix2(N,C,Gp,Ia);

% equilibrium matrix of hinge
[phpn_e,phTpn,theta]=jacobian_ori(node_in_hinge,N,E_n_total);       % jacobian matrix
A_o=[A_2,phTpn];
A_o_a=Ia'*A_o;
% A_o_a=A_o;
%% SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_o_a);         % equilibrium of truss with hinge
% [U1,U2,V1,V2,S1]=tenseg_svd(A_2);           % equilibrium of turss without hinge

%% self-stress design (of truss)
t=zeros(ne,1);      %member force
q=t./l;             % force density
%% self-stress design (of hinge)
M=zeros(n_h,1);
%% cross sectional design (of truss)
A_c=1e-4*ones(ne,1);
A_c(n_b+4*p+1:end)=1e-7*ones(2*p,1);        %reduce the area of vertical string
E_c=1e6*ones(ne,1);
index_b=[1:ne]';              % index of bar in compression
index_s=setdiff(1:size(S,1),index_b);	% index of strings
%% hinge section design  (of hinge)
% k_h=1/12*E_c(index_h).*l(index_h)*thick^3;
% k_h(index_rh_in_h)=1e2*1/12*E_c(index_rh).*l(index_rh)*thick^3;      % increase stiffness of rigid hinge
k_h=1/12*E_c(index_h).*l(1)*thick^3;
%% rest length (of truss), initial angle (of hinge)
l0_c=0.9*l;                     %rest length of truss
theta_0=theta;     % initial angle of hinge (give different value)
%% tangent stiffness matrix of bars
% [Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS(Ia,C,S,q,A_1a,E_c,A_c,l_c);
[Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS2(Ia,C,q,A_2ac,E_c,A_c,l0_c);

%% tangent stiffness matrix of hinge

%ph2px2 is the hessian matrix of theta to nodal coordinate
[ph2pn2_e,ph2pn2]=hessian_ori(node_in_hinge,N,E_n);         % calculate hessian matrix
G=cell2mat(ph2pn2);

%% tangent stiffness of the whole origami

K_t_oa=Kt_aa+Ia'*(phTpn*diag(k_h)*phTpn'+G*kron(M,eye(3*nn)))*Ia;

[K_mode,D1] = eig(K_t_oa);         % eigenvalue of tangent stiffness matrix
k=diag(D1); 
[k, ind] = sort(k);
K_mode = K_mode(:, ind);
% plot the mode shape of tangent stiffness matrix
num_plt=1:9;
plot_mode_ori(round(K_mode,12),k,N,Ia,C_bar,C_s,C_rot_h,C_rig_h,l,'tangent stiffness matrix',...
    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)',num_plt,0.3,saveimg,3,Ca);

%% mass matrix and damping matrix
rho=1;
mass=rho.*A_c.*l0_c;
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0;     %damping coefficient
D=d*2*max(sqrt(mass.*E_c.*A_c./l0_c))*eye(3*nn);    %critical damping
%% vibration mode analysis
[V_mode,D1] = eig(K_t_oa,Ia'*M*Ia);         % calculate vibration mode
w_2=diag(D1);                                    % eigen value of 
omega=real(sqrt(w_2))/2/pi;                   % frequency in Hz
if 0
plot_mode_ori(K_mode,k,N,Ia,C_b,C_s,[],C_rig_h,l,'natrual vibration',...
    'Order of Vibration Mode','Frequency (Hz)',num_plt,0.2,saveimg,3,Ca);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Nonlinear statics analysis %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% external force, forced motion of nodes, shrink of strings
% calculate external force and 
ind_w=[];w=[];   %external force in Z 
ind_dnb=[]; dnb0=[];
ind_dl0_c=[n_b+[1:4*p]]'; dl0_c=-0.5*l0_c(ind_dl0_c);
ind_theta_0=[]'; dtheta_0=[]';        % initial angel change with time
% [w_t,dnb_t,l0_ct,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,l0_c,b,gravity,[0;9.8;0],C,mass);
[w_t,dnb_t,l0_ct,theta_0_t,Ia_new,Ib_new]=tenseg_load_prestress_ori(substep,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,ind_theta_0,dtheta_0,theta_0,l0_c,b,gravity,[0;9.8;0],C,mass);

% % add initial defect by first mode shape
% N_d=N+0.0*max(l)*reshape(Ia*K_mode(:,1),3,[]);
%% Step1: equilibrium calculation
% input data
data.N=N; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia_new; data.Ib=Ib_new;data.S=S;
data.E=E_c; data.A=A_c; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;  % external force
data.dnb_t=dnb_t;           % forced movement of pinned nodes
data.l0_t=l0_ct;            % forced change of rest length
data.theta_0_t=theta_0_t;   % forced change of initial angle
data.k_h=k_h;               % stiffness of hinge
data.E_n=E_n;               % transfer matrix from matrix to structure
data.node_in_hinge=node_in_hinge;       % node in triangle element in hinge
data.substep=substep;    % substep
data.InitialLoadFactor=0.001;
data.MaxIcr=1000;
data.LoadType='Substep'; % 'Force' or 'Displacement' or 'Substep'
data.StopCriterion=@(U)(norm(U)>20);


% nonlinear analysis
data_out1=static_solver_ori_2(data);

Fhis=data_out1.Fhis;          %load factor
t_t=data_out1.t_out;          %member force in every step
n_t=data_out1.n_out;          %nodal coordinate in every step
M_out=data_out1.M_out;
l_out=data_out1.l_out;
icrm=size(n_t,2);               % increment
% N_out=data_out1.N_out;
%% plot member force 
tenseg_plot_result(Fhis,t_t([n_b+1,ne],:),{'horizontal string','vertical string'},{'Load factor','Force / N'},'plot_member_force.png',saveimg);
% %% plot member length 
% tenseg_plot_result(Fhis,l_out,{},{'Load factor','length / m'},'plot_member_length.png',saveimg);

%% plot hinge moment
tenseg_plot_result(Fhis,M_out(index_roth_in_h,:),{'hinge 1','hinge 2','hinge 3','hinge 4','hinge 5'},{'Load factor','Moment / N \times m'},'plot_hinge_moment.png',saveimg);

%% Plot nodal coordinate curve X Y
tenseg_plot_result(Fhis,n_t(3*[(1:p)*2+1]-2,:),{'3X','5X','7X'},{'Load factor','Coordinate /m)'},'plot_coordinate.png',saveimg);
tenseg_plot_result(Fhis,n_t(3*[(1:p)*2],:),{'2Z','4Z','6Z'},{'Load factor','Coordinate /m)'},'plot_coordinate.png',saveimg);

%% Plot final configuration
num_t=4
 j=linspace(1e-5,1,num_t);
for i=1:num_t
    num=ceil(j(i)*size(n_t,2));
tenseg_plot_ori(reshape(n_t(:,num),3,[]),C_b,C_s,[],C_rig_h,[],[],[38,30],[] ,[],Ca);
%  axis off;
end
%% make video of the dynamic
name=['shelter'];
% tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,R3Ddata);
% tenseg_video_slack(n_t,C_b,C_s,l0_ct,index_s,[],[],[],min(substep,50),name,savevideo,material{2})
tenseg_video_ori(n_t,C_b,C_s,[],C_rig_h,Ca,[],min(icrm,50),name,savevideo,[])
%% save output data
if savedata==1
    save (['shelter','.mat']);
end
