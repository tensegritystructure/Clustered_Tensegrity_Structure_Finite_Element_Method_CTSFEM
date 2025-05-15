function data_out=static_solver_RDT(data)
%solve nonlinear equilibrium equations using modified Newton method
%converge to stable equilibrium, considering substep, for Roller-Driven
%Tensegrity

global A E l0 E_qa E_qb C qa qb dqa  w  f_int l_int nn
% minimize total energy? (1: use, 0: not use) it's time consuming
use_energy=0;
tol = 1e-8;MaxIter = 150; 
%% input data
C=data.C;
ne=data.ne;
nn=data.nn;
E_qa=data.E_qa;
E_qb=data.E_qb;
E_na=data.E_na;
E_sa=data.E_sa;
E_nb=data.E_nb;
E_sb=data.E_sb;
% E_qa=blkdiag(E_na,kron(E_sa,eye(2)));
% E_qb=blkdiag(E_nb,kron(E_sb,eye(2)));
% S=data.S;
index_b=data.index_b;
index_s=data.index_s;
substep=data.substep;
E0=data.E;
consti_data=data.consti_data;
material=data.material;
A=data.A;
% w0=data.w;
if  isfield(data,'w_t')
    if size(data.w_t,2)==substep
        w_t=data.w_t;
    elseif size(data.w_t,2)==1
        w_t=data.w_t*linspace(0,1,substep);
    end
else
    w_t=linspace(0,0,substep);
end

% dXb=data.dXb;
if  isfield(data,'dqb_t')
    if size(data.dqb_t,2)==substep
        dqb_t=data.dqb_t;
    elseif size(data.dqb_t,2)==1
        dqb_t=data.dqb_t*linspace(0,1,substep);
    end
else
    dqb_t=linspace(0,0,substep);
end
% l0_0=data.l0;
if size(data.l0_t,2)==substep
    l0_t=data.l0_t;
elseif size(data.l0_t,2)==1
    l0_t=data.l0_t*linspace(1,1,substep);
end

if  isfield(data,'subsubstep')
    subsubstep=data.subsubstep;
else
    subsubstep=30;          %default ssubstep
end

q0=data.q;
data_out=data;     %initialize output data
data_out.E_out=E0*ones(1,substep);


%% calculate equilibrium
qb0=E_qb\q0;           %pinned node
q=q0;               %initialize configuration
E=E0;
% lamda=linspace(0,1,substep);    %coefficient for substep
num_slack=ne*zeros(substep+1,1);    %num of string slack
qa=E_qa\q;
 u=1e-1;


    MaxIcr =substep;                   
%     b_lambda = data.InitialLoadFactor;          
%     Uhis = zeros(3*nn,MaxIcr);        
%     FreeDofs = find(sum(Ia,2));
    lmd = 0; icrm = 0; %MUL = [U,U];
    Fhis = zeros(MaxIcr,1);
    data_out.t_t=zeros(ne,MaxIcr);        %output member force
    data_out.l_t=zeros(ne,MaxIcr);                % member length
    data_out.q_t=zeros(numel(q),MaxIcr);                % generalized coordinate
    data_out.Kt_aa_t=cell(1,MaxIcr);       %tangent stiffness of truss
%     data_out.K_t_oa_out=cell(1,MaxIcr);       %tangent stiffness of whole struct.
   
    F=w_t(:,end);
 while icrm<MaxIcr %&& ~data.StopCriterion(U)&& lmd<=1 
        icrm = icrm+1;
        iter = 0; err = 1;
        lmd=icrm/MaxIcr;
        Fhis(icrm) = lmd;
        w=w_t(:,icrm);               %external force
        qb=qb0+dqb_t(:,icrm);         %forced node displacement
        l0=l0_t(:,icrm);         %forced enlongation of string
%         theta_0=theta_0_t(:,icrm);     % initial angle

        fprintf('\n icrm = %d, lambda = %6.4f\n',icrm,lmd);
         q=[E_qa,E_qb]*[qa;qb];
         sld=q(3*nn+1:end);
         N=reshape(q(1:3*nn),3,[]);
        n_a=pinv(E_na)*N(:);
        n_b=pinv(E_nb)*N(:);
        sld_a=E_sa'*sld;
        sld_b=E_sb'*sld;
        sld_a_jia=0.5.*[abs(sld_a)+sld_a];
        sld_a_jian=0.5.*[abs(sld_a)-sld_a];
        sld_b_jia=0.5.*[abs(sld_b)+sld_b];
        sld_b_jian=0.5.*[abs(sld_b)-sld_b];
        sld_a_s=[sld_a_jia;sld_a_jian];
        sld_b_s=[sld_b_jia;sld_b_jian];
        qa_s=[n_a;sld_a_s];
        while err>tol && iter<MaxIter
            iter = iter+1;
%% equilibrium & tangent stiffness matrix
            sld_a_s=qa_s(end-1:end);
            n_a=qa_s(1:end-2);
            N(:)=E_na*n_a+E_nb*n_b;
            N=reshape(N(:),3,[]);
            H=N*C';
            sld=[E_sa,-E_sa]*sld_a_s+[E_sb,-E_sb]*sld_b_s;
            q=[N(:);sld];
            l=sqrt(sum(H.^2))'; %bar length
            M_u=0.15;%滑轮和绳索之间的摩擦系数
            Theta_1=2*asin(l(1)/(2*l(3)));%滑轮2,3,5的摩擦角
            Theta_2=pi-Theta_1;%滑轮1,4的摩擦角
            Alpha_1=exp(-Theta_1*M_u);
            Alpha_2=exp(-Theta_2*M_u);
            Alpha_all=[Alpha_2;Alpha_1;Alpha_1;Alpha_2;Alpha_1];%mu_i为第i个接触节点的摩擦系数，theta_i为第i个接触节点的总接触角
            a_g=[0 0 9.8];%a_g为[ax ay az]
            rho=7870;
            C_x=max(C,0);
            C_s=C_x-C;

            % member force (of truss)
            %         q=E.*A.*(1./l0-1./l);      %force density
            l0s=l0-C*sld;
            strain=(l-l0s)./l0s;        %strain of member
            [E,sigma]=stress_strain(consti_data,index_b,index_s,strain,material);
            t=sigma.*A;         %member force
            
           
            % equilibrium matrix (of truss)
            Cell_H=mat2cell(H,3,ones(1,size(H,2)));          % transfer matrix H into a cell: Cell_H
            A_2=kron(C',eye(3))*blkdiag(Cell_H{:})*diag(l.^-1);     % equilibrium matrix
            Kn=kron(C'*diag(l.\t)*C,eye(3));       %stiffness matrix of truss
            epsi=ones(2,1)*1e-14;
            M_jia=E_sa'*(pinv(diag(Alpha_all))*C_x'-C_s')*diag(E)*diag(A)*pinv(diag(l0s))*C;
            M_jian=E_sa'*(pinv(diag(Alpha_all))*C_s'-C_x')*diag(E)*diag(A)*pinv(diag(l0s))*C;
            M_a_s=[M_jia;M_jian]*[E_sa -E_sa];
            m_e=rho.*diag(A)*(l0-C*sld);
            m_p=zeros(5,1);
            g_n=kron((0.5.*abs(C)'*m_e+m_p),a_g');
            IFa1=E_na'*Kn*N(:)-(-E_na'*g_n);
            if floor(log10(max(abs(IFa1(:)))))<0
            IFa1_Num=10^(floor(log10(max(abs(IFa1(:)))))+1);
            else IFa1_Num=10^floor(log10(max(abs(IFa1(:)))));end
            g_s=0.5*rho.*kron((C'*diag(A)*abs(C)),a_g)*N(:);
            sld_a=E_sa'*sld;
            sld_b=E_sb'*sld;
            sld_a_jia=0.5.*[abs(sld_a)+sld_a];
            sld_a_jian=0.5.*[abs(sld_a)-sld_a];
            sld_b_jia=0.5.*[abs(sld_b)+sld_b];
            sld_b_jian=0.5.*[abs(sld_b)-sld_b];
            sld_a_s=[sld_a_jia;sld_a_jian];
            sld_b_s=[sld_b_jia;sld_b_jian];
            zeta_sa_s=[E_sa'*(pinv(diag(Alpha_all))*C_x'-C_s')*t-(-E_sa'*g_s);...
                E_sa'*(pinv(diag(Alpha_all))*C_s'-C_x')*t+(-E_sa'*g_s)];
%             if zeta_sa_s(1)>0
%                 dqa_s(5)=0;end
%             if zeta_sa_s(2)>0
%                 dqa_s(6)=0;end
            X1=max(abs(zeta_sa_s(:)))*ones(2,1)/IFa1_Num;
            if max(abs(sld_a_s(:)))==0
                X2=ones(2,1);
            else    X2=1/(max(abs(sld_a_s(:))))*ones(2,1)*IFa1_Num;end
            PSI_s=((X1.\zeta_sa_s).^2+(X2.*sld_a_s).^2+2*epsi.^2).^0.5-(X1.\zeta_sa_s+X2.*sld_a_s);
            dHs_dqa11=E_na'*(kron((C'*pinv(diag(l))*diag(t)*C),eye(3))+A_2*diag(E)*diag(A)*pinv(diag(l))*A_2')*E_na;
            dHs_dqa12=E_na'*A_2*diag(E)*diag(A)*pinv(diag(l0s))*C*E_sa*[eye(size(E_sa,2)),-eye(size(E_sa,2))];%%缺
            dHs_dqa21=(zeta_sa_s./(X1.*((X1.\zeta_sa_s).^2+(X2.*sld_a_s).^2+2*epsi.^2).^0.5)-1./X1).*kron(eye(2),E_sa')*[pinv(diag(Alpha_all))*C_x'-C_s';pinv(diag(Alpha_all))*C_s'-C_x']...
                        *diag(E)*diag(A)*pinv(diag(l0s))*A_2'*E_na;
            dHs_dqa22=(zeta_sa_s/(X1.*((X1.\zeta_sa_s).^2+(X2.*sld_a_s).^2+2*epsi.^2).^0.5)-1./X1)*M_a_s+(X2.*sld_a_s/(((X1.\zeta_sa_s).^2+(X2.*sld_a_s).^2+2*epsi.^2).^0.5)-X2);

            
            

            % unbalanced force
            IFa=-[E_na'*Kn*N(:)-(-E_na'*g_n);PSI_s];

%             IF=w-(Ki*q+phTpn*M);                   %unbalanced force
%             IFa=Ia'*(w-(Ki*q+phTpn*M));            %unbalanced force
%             Fp_a=Ia'*Fp;                   %see the norm of unbalanced force
            % tangent stiffness matrix (of truss)
            K_T=[Kn+A_2*diag(E.*A./l0)*A_2'-A_2*diag(l.\t)*A_2',A_2*diag(E.*A./l0s)*C;...
            C'*diag(E.*A./l0s)*A_2',C'*diag(E.*A./l0s)*C];
            K_T=0.5*(K_T+K_T');
            K_Taa=E_qa'*K_T*E_qa;
            K_Tab=E_qa'*K_T*E_qb;
            Km=[dHs_dqa11,dHs_dqa12;dHs_dqa21,dHs_dqa22];
            
        if 0
        [V_mode,D]=eig (Km);                       %刚度矩阵特征根
        d=diag(D);                            %eigen value
        lmd=min(d);                     %刚度矩阵最小特征根
        if lmd>0
            Km= Km+u*eye(size( Km)); %修正的刚度矩阵
        else
            Km= Km+(abs(lmd)+u)*eye(size( Km));
        end
        end
           
            dqa_s= Km\IFa;

        x=1;
        % line search
        if use_energy==1
            opt=optimset('TolX',1e-5);
            [x,V]=fminbnd(@energy_RDT,0,1e1,opt);
        end
        


            err = norm(IFa);
            qa_s=qa_s+x*dqa_s;
            fprintf(' iter = %d, err = %6.4e\n',iter,err);
            if err > 1e12, disp('Divergence!'); break; end
        end

    data_out.q_t(:,icrm)=q;
    data_out.l_t(:,icrm)=l;
    data_out.t_t(:,icrm)=t;      %member force
data_out.Kt_t{:,icrm}=K_T;     % tangent stiffness matrix
 end