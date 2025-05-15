clc;
clear;
close all;

%%
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
substep=1;                                     %荷载子步
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=1;             % save data or not (1) yes (0)no
savevideo=1;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
% move_ground=0;          % for earthquake, use pinned nodes motion(1) or add inertia force in free node(0) 
%% %% N C of the structure
% rate1=19.319/50;7
rate1=0.2;       %内环半径与外环半径之比
R=0.8;          %radius
r0=rate1*R;     %内环半径
fd1 = 200;
p=6;          %complexity for cable dome(Outer ring node)；2的倍数
gr_num=1;      %同一个环，组数划分（1，2，3，4，6）；节点数需为其的整数倍数

    h3=-0.4;       %下环高度
    h2=0.4;       %上环高度
    h4=h3*r0/cos(pi/6)/R;         %下中环高度
    h1=r0*(h2-h4)/((h2-h4)^2+(R*cos(pi/6)-r0)^2)^0.5+h4;       %上中环高度    
  
    % generate node in one unit
    beta1=pi/(p); beta2=4*pi/p;         %  两点间旋转量
    beta3=2*pi/p;         %  整体角度
    T1=[cos(beta1) -sin(beta1) 0
    sin(beta1) cos(beta1) 0
    0 0 1];            %内节点旋转量

    T2=[cos(beta2) -sin(beta2) 0
    sin(beta2) cos(beta2) 0
    0 0 1];            %上下节点旋转量

    T3=[cos(beta3) -sin(beta3) 0
    sin(beta3) cos(beta3) 0
    0 0 1];


    N_1_0=[T1*[r0+r0*(R*cos(pi/6)-r0)/((h2-h4)^2+(R*cos(pi/6)-r0)^2)^0.5;0;h1]];       %上内初节点
    N_4_0=[T2*[r0/cos(pi/6);0;h4]];       %下内初节点
 
    N_2_0=[T2*[R;0;h2]];      %上初节点
    N_3_0=[T2*[R;0;h3]];      %下初节点

    N_1=[];
    N_2=[];
    N_3=[];
 N_4=[];
 
    for i=1:p  %内环节点
     N_1=[N_1,T3^(i-1)*N_1_0];
    end
    for i=0:1:p  %内环节点
        N_4=[N_4,T3^(i-1)*N_4_0];
    end
    for i=0:p-1    %上环节点
     N_2=[N_2,T3^(i-1)*N_2_0];
    end
    for i=0:p-1    %下环节点
     N_3=[N_3,T3^(i-1)*N_3_0];
    end

    N=[N_3,N_1,N_2,N_4];  

    C_b_in=[];

    C_s_in=[[1:1:p]',[3*p+1:1:4*p]';
        [2*p+1:1:3*p]',[p+1:1:2*p]';[2*p+1:1:3*p]',[p+2:1:2*p,p+1]';
        [p+1:1:2*p]',[4*p,3*p+1:1:4*p-1]';[p+1:1:2*p]',[3*p+1:1:4*p]';
         [3*p+1:1:4*p]',[3*p+2:1:4*p,3*p+1]'];

    C_b = tenseg_ind2C(C_b_in,N);  
    C_s = tenseg_ind2C(C_s_in,N);   
    C=[C_s;C_b];
    [ne,nn]=size(C);        % ne:No.of element;nn:No.of node
    % Plot the structure to make sure it looks right
    
    tenseg_plot(N,C_b,C_s);
    hold on;
    title('Cable dome');
view(30,0)
%% %% Boundary constraints
% pinned_X=([1:1:2*p,2*p+1,3*p+1])'; pinned_Y=([1:1:2*p,2*p+1,3*p+1])'; pinned_Z=([1:1:2*p,2*p+1,3*p+1])';

pinned_X=([1:1:p,2*p+1:1:3*p])'; pinned_Y=([1:1:p,2*p+1:1:3*p])'; pinned_Z=([1:1:p,2*p+1:1:3*p])';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);
%% %% Group/Clustered

% information 
%generate group index
% gr=[];
[gr] = cable_ring_gr(gr_num,p);
Gp=tenseg_str_gp3(gr,C);   %generate group matrix
% S=eye(ne);                  % no clustering matrix

S=Gp';                      % clustering matrix is group matrix

tenseg_plot_CTS(N,C,[],S);
%% %% self-stress design
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N,C,Gp,Ia);
A_1ac=A_1a*S';          %equilibrium matrix CTS
A_2ac=A_2a*S';          %equilibrium matrix CTS
l_c=S*l;                % length vector CTS
%SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_1ag);

%external force in equilibrium design
w0=zeros(numel(N),1); w0a=Ia'*w0;

%prestress design
%index_gp=[1,3]/[2]; % number of groups with designed force

switch gr_num
    case 1
index_gp=[3]; %一组
fd=[fd1];
    case 2
index_gp=[1,2]; %二组
fd=fd1*ones(2,1);
    case 3
index_gp=[1,2,3]; %三组
fd=fd1*ones(3,1);
    case 4
index_gp=[5,6,7,8]; %四组
fd=fd1*ones(4,1);
    case 6
index_gp=[13,14,15,16,17,18];%六组
fd=fd1*ones(6,1);
    case 12
index_gp=[4*p+1:5*p];%六组
fd=fd1*ones(6,1);
end 

%fd=[1000; 1000];                       % force in bar is given as -1000
% fd=[1000];
[q_gp,t_gp,q,t]=tenseg_prestress_design(Gp,l,l_gp,A_1ag,V2,w0a,index_gp,fd);    %prestress design
t_c=pinv(S')*t;
q_c=pinv(S')*q;
%% cross sectional design
index_b=find(t_c<0);              % index of bar in compression
index_s=setdiff(1:size(S,1),index_b);	% index of strings

[A_b,A_s,A_c,A,r_b,r_s,r_gp,radius,E_c,l0_c,rho,mass_c]=tenseg_minimass(t_c,l_c,eye(size(S,1)),sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,thick,hollow_solid);
E=S'*E_c;     %Young's modulus CTS
 
A=S'*A_c;     % Cross sectional area CTS,所有成员的截面积
l0=(t+E.*A).\E.*A.*l;
mass=S'*rho.*A.*l0;
% % Plot the structure with radius
% R3Ddata.Bradius=interp1([min(radius),max(radius)],[0.03,.1],r_b);
% R3Ddata.Sradius=interp1([min(radius),max(radius)],[0.03,.1],r_s);
% R3Ddata.Nradius=0.1*ones(nn,1);
% tenseg_plot(N,C_b,C_s,[],[],[],'Double layer prism',R3Ddata);
%% tangent stiffness matrix
% ft_aa=20*ones(36,36);
[Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS(Ia,C,S,q,A_1a,E_c,A_c,l_c);
% plot the mode shape of tangent stiffness matrix
% ut_aa=ft_aa/Kt_aa
% [u_mode,D2] = eig(ut_aa);  
% u=diag(D2);   
lmd_k = min(k);
num_plt=1:4;
if 1
plot_mode(K_mode,k,N,Ia,C_b,C_s,l,'tangent stiffness matrix',...
    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)',num_plt,0.8,saveimg,3);
end
%% input file of ANSYS
% ansys_input_gp(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),'tower');

%% mass matrix and damping matrix
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0;     %damping coefficient
D=d*2*max(sqrt(mass.*E.*A./l0))*eye(3*nn);    %critical damping
%% mode analysis
[V_mode,D1] = eig(Kt_aa,Ia'*M*Ia);         % calculate vibration mode
w_2=diag(D1);                                    % eigen value of 
omega=real(sqrt(w_2))/2/pi;                   % frequency in Hz
if 0
plot_mode(V_mode,omega,N,Ia,C_b,C_s,l,'natrual vibration',...
    'Order of Vibration Mode','Frequency (Hz)',num_plt,0.8,saveimg,3);
end

%% external force, forced motion of nodes, shrink of strings
% calculate external force and 
ind_w=[];w=[];
ind_dnb=[]; dnb0=[];
ind_dl0_c=[3]; dl0_c=[-0.5];
[w_t,dnb_t,l0_ct,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,l0_c,b,gravity,[0;0;9.8],C,mass);
% [w_t,dnb_t,l0_t,Ia,Ib]=tenseg_load_prestress(substep,ind_w,w0,ind_dn,dn0,ind_l0,dl0,l0,b,gravity,acc,C,mass)

%% Step1: equilibrium calculation
% input data
data.N=N; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia_new; data.Ib=Ib_new;data.S=S;
data.E=E_c; data.A=A_c; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;  % external force
data.dnb_t=dnb_t;% forced movement of pinned nodes
data.l0_t=l0_ct;% forced movement of pinned nodes
data.substep=substep;    % substep

% nonlinear analysis
% data_out=static_solver(data);        %solve equilibrium using mNewton method
% data_out=static_solver2(data);        %solve equilibrium using mNewton method
data_out1=static_solver_CTS(data);
% data_out{i}=equilibrium_solver_pinv(data);        %solve equilibrium using mNewton method

t_t=data_out1.t_out;          %member force in every step  (每根绳索的力）
n_t=data_out1.n_out;          %nodal coordinate in every step
lmd_out=data_out1.lmd_out;  %最小刚度特征值
N_out=data_out1.N_out;
tenseg_plot( N_out{:},C_b,C_s,[],[],[])
   % tenseg_plot( N,C_b,C_s,fig_handle,highlight_nodes,view_vec, PlotTitle, R3Ddata)
tenseg_plot_CTS(N_out{:},C,[],S);
%%   stiffness 
N2=N_out{end};         % final configuration 
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N2,C,Gp,Ia);
A_1ac=A_1a*S';          % equilibrium matrix CTS
A_2ac=A_2a*S';          % equilibrium matrix CTS
l_c=S*l;                % length vector CTS

t=t_t(:,end);           % members' force
q=t./l;                 % members' force density
%% tangent stiffness matrix analysis
[Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS(Ia,C,S,q,A_1a,E_c,A_c,l_c);
% plot the mode shape of tangent stiffness matrix
num_plt=1:4;
if 1
plot_mode(K_mode,k,N2,Ia,C_b,C_s,l,'tangent stiffness matrix',...
    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)',num_plt,0.8,saveimg,[25,30]);
end
k1=min(k);
omega1=min(omega);
%% mode analysis
[V_mode,D1] = eig(Kt_aa,Ia'*M*Ia);         % calculate vibration mode
w_2=diag(D1);                                    % eigen value of 
omega=real(sqrt(w_2))/2/pi;                   % frequency in Hz
if 0
plot_mode(V_mode,omega,N2,Ia,C_b,C_s,l,'natrual vibration',...
    'Order of Vibration Mode','Frequency (Hz)',num_plt,0.8,saveimg,[25,30]);
end
%% Step 2: change rest length of strings

substep=10;
ind_dnb=[]; dnb0=[];

dl0_c_1 = -250;  %内环收缩量
dl0_c_2 = 200;   %上环伸长量
dl0_c_3 = 200;   %下环伸长量
% 
[ind_dl0_c,dl0_c] = tenseg_dl0_1(gr_num,dl0_c_1,dl0_c_2,dl0_c_3)

% ind_dl0_c=tenseg_dl0(gr_num);


[w_t,dnb_t,l0_ct,Ia_new1,Ib_new1]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,l0_c,b,gravity,[0;9.8;0],C,mass);
       % [w_t,dnb_t,l0_t,Ia,Ib]=tenseg_load_prestress(substep,ind_w,w0,ind_dn,dn0,ind_l0,dl0,l0,b,gravity,acc,C,mass)
data.w_t=w_t;  % external force
data.dnb_t=dnb_t;% forced movement of pinned nodes
data.l0_t=l0_ct;% forced movement of pinned nodes
data.N=N_out{end};
data.substep=substep;    % substep

data_out=static_solver_CTS(data);

t_t=data_out.t_out;          %member force in every stepn_t=data_out.n_out;          %nodal coordinate in every step
N_out=data_out.N_out;
K_out=data_out.K_out;
l_c_out=data_out.l_c_out;   %各步态下时的弦长长度

%% plot member force 
tenseg_plot_result(1:substep,t_t([3*p,5*p,p],:),{'上斜索','内斜索','下环索'},{'Load step','Force (N)'},'plot_member_force.png',saveimg);


%% Plot nodal coordinate curve X Y
tenseg_plot_result(1:substep,n_t([3*2-2,3*2],:),{'2X','2Z'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);

%% Plot final configuration
% tenseg_plot_catenary( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,0],[],[],l0_ct(index_s,end))
% tenseg_plot( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[])
tenseg_plot( reshape(n_t(:,10),3,[]),C_b,C_s,[],[],[])

%% save output data
if savedata==0
    save (['cable_ring',material{1},'.mat']);
end
%% make video of the dynamic
name=['cable_ring'];
% tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,R3Ddata);
% tenseg_video_slack(n_t,C_b,C_s,l0_ct,index_s,[],[],[],min(substep,50),name,savevideo,material{2})
tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,material{2})

%output data to tecplot
%tenseg_tecplot(C,n_t,t_t,interp1([min(radius),max(radius)],[0.2,0.8],radius));