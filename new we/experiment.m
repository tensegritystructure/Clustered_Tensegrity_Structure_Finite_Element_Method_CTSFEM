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
%  rate1=0.3;       %内环半径与外环半径之比
xp=1;                %细分开合度用
 rate0=linspace(0.3,0.3,xp);
p=6;          %complexity for cable dome(Outer ring node)；2的倍数
n_t0=zeros(3*11*p,xp);%7*p为点数
l0_cxp=zeros(11*p,xp);%10*p为索和杆的总数
t_cxp=zeros(10*p,xp);
t_xp=zeros(16*p,xp);
A_xp=zeros(16*p,xp);
l_xp=zeros(16*p,xp);
l_cxp=zeros(10*p,xp);
l_0xp=zeros(11*p,xp);
q_xp=zeros(16*p,xp);
A_1axp=zeros(16*p*p*18,xp);
for ii=1:xp
    rate1=rate0(:,ii);
R=0.27;          %radius
r0=rate1*R;     %内环半径
fd1 = 200;
p=6;          %complexity for cable dome(Outer ring node)；2的倍数
gr_num=6;      %同一个环，组数划分（1，2，3，4，6）；节点数需为其的整数倍数
sgg=0.002;      %三角形杆的高
gbj=(2*pi-pi*(p-2)/p-pi/3)/2; %小三角杆旁边的角的角度

    h3=-0.14;       %下环高度
    h2=0.14;       %上环高度
    h4=h3*r0/R;         %下中环高度
    h1=(((r0*cos(pi/p))^2+h4^2)^0.5)*(h2-h4)/((h2-h4)^2+(R*cos(pi/p)-r0*cos(pi/p)-sgg/cos(pi/6)*sin(gbj))^2)^0.5+h4;       %上中环高度    
    h5=(sgg/cos(pi/6))*(h2-h4)/((h2-h4)^2+(R*(cos(pi/p))^2-r0)^2)^0.5+h4;       %上中环高度
    h6=(sgg/cos(pi/6))*(h3-h4)/((h3-h4)^2+(R*(cos(pi/p))^2-r0)^2)^0.5+h4;
    
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


    N_1_0=[T1*[r0*cos(pi/p)+sgg/cos(pi/6)*sin(gbj)+(R*cos(pi/p)-r0*cos(pi/p)-sgg/cos(pi/6)*sin(gbj))*(h1-h4)/(h2-h4);0;h1]];       %上内初节点
    N_5_0=[T2*[r0;0;h4]];       %下内初节点
    N_6_0=[T2*[r0+sgg;sgg*tan(pi/6);h4]];       %下内外初节点1
    N_7_0=[T2*[r0+sgg;-sgg*tan(pi/6);h4]];       %下内外初节点2
    N_2_0=[T2*[R;0;h2]];      %上初节点
    N_3_0=[T2*[R;sgg*tan(pi/6);h3]];      %下初节点1
    N_4_0=[T2*[R;-sgg*tan(pi/6);h3]];      %下初节点1
    N_8_0=[T1*[R*cos(pi/p);0;h2]];    %环索主力上斜索
    N_9_0=[T1*[R*cos(pi/p);0;h3]];   %环索主力下斜索
    N_10_0=[T2*[(sgg/cos(pi/6))*(R*(cos(pi/p))^2-r0)/((h2-h4)^2+(R*(cos(pi/p))^2-r0)^2)^0.5+r0;0;h5]];   %环索主力上斜索杆
    N_11_0=[T2*[(sgg/cos(pi/6))*(R*(cos(pi/p))^2-r0)/((h3-h4)^2+(R*(cos(pi/p))^2-r0)^2)^0.5+r0;0;h6]];   %环索主力下斜索杆

    N_1=[];
    N_2=[];
    N_3=[];
    N_4=[];
    N_5=[];
    N_6=[];
    N_7=[];
    N_8=[];
 	N_9=[];
    N_10=[];
    N_11=[];
 
    for i=1:p  %上内初节点
     N_1=[N_1,T3^(i-1)*N_1_0];
    end
    for i=0:p-1  %下内初节点
        N_5=[N_5,T3^(i-1)*N_5_0];
    end
     for i=0:p-1  %下内外初节点1
        N_6=[N_6,T3^(i-1)*N_6_0];
     end
    for i=0:p-1  %下内外初节点2
        N_7=[N_7,T3^(i-1)*N_7_0];
    end
    for i=0:p-1    %上环节点
     N_2=[N_2,T3^(i-1)*N_2_0];
    end
    for i=0:p-1    %下环节点1
     N_3=[N_3,T3^(i-1)*N_3_0];
    end
    for i=0:p-1    %下环节点2
     N_4=[N_4,T3^(i-1)*N_4_0];
    end
    for i=0:p-1    %下环节点2
     N_8=[N_8,T3^(i-1)*N_8_0];
    end
    for i=0:p-1    %下环节点2
     N_9=[N_9,T3^(i-1)*N_9_0];
    end
    for i=0:p-1    %下环节点2
     N_10=[N_10,T3^(i-1)*N_10_0];
    end
    for i=0:p-1    %下环节点2
     N_11=[N_11,T3^(i-1)*N_11_0];
    end
    N=[N_3,N_4,N_1,N_2,N_5,N_6,N_7,N_8,N_9,N_10,N_11];  

    C_b_in=[[4*p+1:1:5*p]',[5*p+1:1:6*p]';[4*p+1:1:5*p]',[6*p+1:1:7*p]';[5*p+1:1:6*p]',[6*p+1:1:7*p]';...
        [9*p+1:1:10*p]',[4*p+1:1:5*p]';[10*p+1:1:11*p]',[4*p+1:1:5*p]'];

    C_s_in=[[1:1:p]',[5*p+1:1:6*p]';[p+1:1:2*p]',[6*p+1:1:7*p]';[5*p+1:1:6*p]',[2*p+2:1:3*p,2*p+1]';...
        [6*p+1:1:7*p]',[2*p+1:1:3*p]';[2*p+1:1:3*p]',[3*p+1:1:4*p]';[2*p+2:1:3*p,2*p+1]',[3*p+1:1:4*p]';[4*p+1:1:5*p]',[4*p+2:1:5*p,4*p+1]';...
        [7*p+2:1:8*p,7*p+1]',[9*p+1:1:10*p]';[7*p+3:1:8*p,7*p+1,7*p+2]',[9*p+1:1:10*p]';[8*p+2:1:9*p,8*p+1]',[10*p+1:1:11*p]';[8*p+3:1:9*p,8*p+1,8*p+2]',[10*p+1:1:11*p]'];

    C_b = tenseg_ind2C(C_b_in,N);  
    C_s = tenseg_ind2C(C_s_in,N);   
    C=[C_s;C_b];
    Ca=[[2*p+1:3*p;6*p+1:7*p;6*p,5*p+1:6*p-1],[6*p+1:7*p;6*p,5*p+1:6*p-1;p,1:p-1],[6*p+1:7*p;p,1:p-1;p+1:2*p]];
    C_h_in=[[1:1:p]',[5*p+1:1:6*p]';[p+1:1:2*p]',[6*p+1:1:7*p]';[5*p+1:1:6*p]',[2*p+2:1:3*p,2*p+1]';[6*p+1:1:7*p]',[2*p+1:1:3*p]'];
    C_rh_in=[[1:1:p]',[6*p+2:7*p,6*p+1]'];
    C_h = tenseg_ind2C(C_h_in,N);
    C_rh= tenseg_ind2C(C_rh_in,N);
    [ne,nn]=size(C);        % ne:No.of element;nn:No.of node
    % Plot the structure to make sure it looks right
    n_t01=N(:);
    n_t0(:,ii)=n_t01;

    tenseg_plot(N,C_b,C_s);
    hold on;
    title('Cable dome');
view(30,0); 
%% %% Boundary constraints
% pinned_X=([1:1:2*p,2*p+1,3*p+1])'; pinned_Y=([1:1:2*p,2*p+1,3*p+1])'; pinned_Z=([1:1:2*p,2*p+1,3*p+1])';

pinned_X=([1:1:2*p,3*p+1:1:4*p,8*p+1:1:9*p,7*p+1:1:8*p])'; pinned_Y=([1:1:2*p,3*p+1:1:4*p,8*p+1:1:9*p,7*p+1:1:8*p])'; pinned_Z=([1:1:2*p,3*p+1:1:4*p,8*p+1:1:9*p,7*p+1:1:8*p])';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);

%% %% Group/Clustered

%information 
%generate group index
% gr=[];

%   matrix_1 = [5*p+1:8*p];
%   Gr = num2cell(matrix_1)';
% gr={[2*p+1:4*p]';[1:2*p]';[4*p+1:5*p]'};  % 上，下，环；一组
% gr_1 = [gr;Gr];

[gr] = cable_ring_gr_copy(gr_num,p);
Gp=tenseg_str_gp3(gr,C);   %generate group matrix
% S=eye(ne);                  % no clustering matrix

S=Gp';                      % clustering matrix is group matrix

% tenseg_plot_CTS(N,C,[],S);
%% %% self-stress design
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N,C,Gp,Ia);%85-87，16-20
A_1ac=A_1a*S';          %equilibrium matrix CTS，85-87
A_2ac=A_2a*S';          %equilibrium matrix CTS
l_c=S*l;                % length vector CTS，19
%SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_1ag);%U1列空间；V1行空间；U2左零空间；V2零空间；S特征值，(2.49)

%external force in equilibrium design
w0=zeros(numel(N),1); w0a=Ia'*w0; %numel算矩阵的元素数

%prestress design
%index_gp=[1,3]/[2]; % number of groups with designed force

switch gr_num
    case 1
index_gp=[4]; %一组
fd=[fd1];
    case 2
index_gp=[1,2]; %二组
fd=fd1*ones(2,1);
    case 3
index_gp=[1,2,3]; %三组
fd=fd1*ones(3,1);
%     case 4
% index_gp=[5,6,7,8]; %四组
% fd=fd1*ones(4,1);
    case 6
index_gp=[4*p+1];%六组
fd=fd1*ones(1,1);
    case 12
index_gp=[3*p+1:4*p];%六组
fd=fd1*ones(6,1);
end

%fd=[1000; 1000];                       % force in bar is given as -1000
% fd=[1000];
[q_gp,t_gp,q,t]=tenseg_prestress_design(Gp,l,l_gp,A_1ag,V2,w0a,index_gp,fd);    %prestress design
t_c=pinv(S')*t;  %31
q_c=pinv(S')*q;  %33
q_xp(:,ii)=q;
%% cross sectional design
index_b=find(t_c<0);              % index of bar in compression
index_s=setdiff(1:size(S,1),index_b);	% index of strings

[A_b,A_s,A_c,A,r_b,r_s,r_gp,radius,E_c,l0_c,rho,mass_c]=tenseg_minimass(t_c,l_c,eye(size(S,1)),sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,thick,hollow_solid);
l_xp(:,ii)=l;
l_cxp(:,ii)=l_c;
A=S'*A_c;
A_xp(:,ii)=A;
t_cxp(:,ii)=t_c;
t_xp(:,ii)=t;
A_1az=reshape(A_1a,[],1);
A_1axp(:,ii)=A_1az;
end
ratio_max=ratio_max(h2,R,p,sgg,gbj,h3);
N=reshape(n_t0(:,1),3,[]);
A_row=max(A_xp,[],2);
E=S'*E_c;     %Young's modulus CTS 28
t_1=min(t_cxp,[],2);
A_c_row=pinv(S')*A_row;     % Cross sectional area CTS,所有成员的截面积 27
t_cr=(A_row.^2).*pi.*E./(4.*l.^2);%压杆稳定

l0_xp=(t_xp+repmat(E.*A_row,1,xp)).\repmat(E.*A_row,1,xp).*l_xp; %？30+
l0=l0_xp(:,xp);
l0_cxp=S*l0_xp;
l0_c=l0_cxp(:,xp);
mass=S'*rho.*A_row.*l0;   %30 ?
% % Plot the structure with radius
% R3Ddata.Bradius=interp1([min(radius),max(radius)],[0.03,.1],r_b);
% R3Ddata.Sradius=interp1([min(radius),max(radius)],[0.03,.1],r_s);
% R3Ddata.Nradius=0.1*ones(nn,1);
% tenseg_plot(N,C_b,C_s,[],[],[],'Double layer prism',R3Ddata);
tenseg_plot_CTS(N,C,[],S);
%% patch
%绘图
for i=1:1:p-1
patch(N(1,[i,i+p+1,i+6*p+1,i+5*p]),N(2,[i,i+p+1,i+6*p+1,i+5*p]),N(3,[i,i+p+1,i+6*p+1,i+5*p]),2);
patch(N(1,[i+2*p+1,i+5*p,i+6*p+1]),N(2,[i+2*p+1,i+5*p,i+6*p+1]),N(3,[i+2*p+1,i+5*p,i+6*p+1]),2);
end
patch(N(1,[1+2*p,6*p,1+6*p]),N(2,[1+2*p,6*p,1+6*p]),N(3,[1+2*p,6*p,1+6*p]),2);
patch(N(1,[p,p+1,6*p+1,6*p]),N(2,[p,p+1,6*p+1,6*p]),N(3,[p,p+1,6*p+1,6*p]),2);

%% tangent stiffness matrix
% ft_aa=20*ones(36,36);
lmd_kxp=zeros(1,xp);
for jj=1:xp
    q=q_xp(:,jj);
A_1az=A_1axp(:,jj);
A_1a=reshape(A_1az,[],96);
l_c=l_cxp(:,jj);
[Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS(Ia,C,S,q,A_1a,E_c,A_c_row,l_c);
% % plot the mode shape of tangent stiffness matrix
% % ut_aa=ft_aa/Kt_aa
% % [u_mode,D2] = eig(ut_aa);  
% % u=diag(D2);   
lmd_k = min(k);
lmd_kxp(:,jj)=lmd_k;
end
num_plt=1:4;
if 1
plot_mode(K_mode,k,N,Ia,C_b,C_s,l,'tangent stiffness matrix',...
    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)',num_plt,0.8,saveimg,3);
end
% %% input file of ANSYS
% % ansys_input_gp(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),'tower');
% 
% %% mass matrix and damping matrix
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
end    %5月14日
% 
%% external force, forced motion of nodes, shrink of strings
% calculate external force and 
ind_w=[];w=[];
ind_dnb=[]; dnb0=[];
ind_dl0_c=[]; dl0_c=[];
[w_t,dnb_t,l0_ct,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,l0_c,b,gravity,[0;0;9.8],C,mass);
% [w_t,dnb_t,l0_t,Ia,Ib]=tenseg_load_prestress(substep,ind_w,w0,ind_dn,dn0,ind_l0,dl0,l0,b,gravity,acc,C,mass)

%% Step1: equilibrium calculation
% input data
data.N=N; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia_new; data.Ib=Ib_new;data.S=S;
data.E=E_c; data.A=A_c_row; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;  % external force
data.dnb_t=dnb_t;% forced movement of pinned nodes
data.l0_t=l0_cxp;% forced movement of pinned nodes
data.substep=substep;    % substep

% % nonlinear analysis
% % data_out=static_solver(data);        %solve equilibrium using mNewton method
% % data_out=static_solver2(data);        %solve equilibrium using mNewton method
data_out1=static_solver_CTS(data);
% % data_out{i}=equilibrium_solver_pinv(data);        %solve equilibrium using mNewton method

t_t=data_out1.t_out;          %member force in every step  (每根绳索的力）
n_t=data_out1.n_out;%nodal coordinate in every step
% 
% 
% 
% lmd_out=data_out1.lmd_out;  %最小刚度特征值
% N_out=data_out1.N_out;
% tenseg_plot( N_out{:},C_b,C_s,[],[],[]);
%    % tenseg_plot( N,C_b,C_s,fig_handle,highlight_nodes,view_vec, PlotTitle, R3Ddata)
% tenseg_plot_CTS(N_out{:},C,[],S);

% %% Step 2: change rest length of strings
% 
% substep=10;
% ind_dnb=[]; dnb0=[];
% 
% dl0_c_1 = -250;  %内环收缩量
% dl0_c_2 = 200;   %上环伸长量
% dl0_c_3 = 200;   %下环伸长量
% % 
% [ind_dl0_c,dl0_c] = tenseg_dl0_1(gr_num,dl0_c_1,dl0_c_2,dl0_c_3)
% 
% % ind_dl0_c=tenseg_dl0(gr_num);
% 
% 
% [w_t,dnb_t,l0_ct,Ia_new1,Ib_new1]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,l0_c,b,gravity,[0;9.8;0],C,mass);
%        % [w_t,dnb_t,l0_t,Ia,Ib]=tenseg_load_prestress(substep,ind_w,w0,ind_dn,dn0,ind_l0,dl0,l0,b,gravity,acc,C,mass)
% data.w_t=w_t;  % external force
% data.dnb_t=dnb_t;% forced movement of pinned nodes
% data.l0_t=l0_ct;% forced movement of pinned nodes
% data.N=N_out{end};
% data.substep=substep;    % substep
% 
% data_out=static_solver_CTS(data);
% 
% t_t=data_out.t_out;          %member force in every stepn_t=data_out.n_out;          %nodal coordinate in every step
% N_out=data_out.N_out;
% K_out=data_out.K_out;
% l_c_out=data_out.l_c_out;   %各步态下时的弦长长度
% 
% %% plot member force 
% tenseg_plot_result(1:substep,t_t([3*p,5*p,p],:),{'上斜索','内斜索','下环索'},{'Load step','Force (N)'},'plot_member_force.png',saveimg);
% 
% 
% %% Plot nodal coordinate curve X Y
% tenseg_plot_result(1:substep,n_t([3*2-2,3*2],:),{'2X','2Z'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
% 
% %% Plot final configuration
% % tenseg_plot_catenary( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,0],[],[],l0_ct(index_s,end))
% % tenseg_plot( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[])
% tenseg_plot( reshape(n_t(:,10),3,[]),C_b,C_s,[],[],[])
% 
% %% save output data
% if savedata==0
%     save (['cable_ring',material{1},'.mat']);
% end
%% make video of the dynamic
name=['cable_ring1'];
% tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,R3Ddata);
% tenseg_video_slack(n_t,C_b,C_s,l0_ct,index_s,[],[],[],min(substep,50),name,savevideo,material{2})
%tenseg_video(n_t0,C_b,C_s,[],xp,name,savevideo,material{2});
tenseg_video_ori(n_t0,C_b,C_s,C_h,C_rh,Ca,[],xp,name,savevideo,material{2})
%output data to tecplot
%tenseg_tecplot(C,n_t,t_t,interp1([min(radius),max(radius)],[0.2,0.8],radius));
