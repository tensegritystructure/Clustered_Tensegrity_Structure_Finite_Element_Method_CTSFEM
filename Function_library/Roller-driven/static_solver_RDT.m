function data_out=static_solver_RDT(data)
%solve nonlinear equilibrium equations using modified Newton method
%converge to stable equilibrium, considering substep, for Roller-Driven
%Tensegrity

global A E l0 E_qa E_qb C qa qb dqa  w  f_int l_int
% minimize total energy? (1: use, 0: not use) it's time consuming
use_energy=0;
tol = 1e-6;
%% input data
C=data.C;
ne=data.ne;
nn=data.nn;
E_qa=data.E_qa;
E_qb=data.E_qb;
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
        dXb_t=data.dqb_t;
    elseif size(data.dqb_t,2)==1
        dXb_t=data.dqb_t*linspace(0,1,substep);
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
qb0=E_qb'*q0;           %pinned node
q=q0;               %initialize configuration
E=E0;
% lamda=linspace(0,1,substep);    %coefficient for substep
num_slack=ne*zeros(substep+1,1);    %num of string slack
qa=E_qa'*q;
 u=1e-1;


    MaxIcr =substep;                   
%     b_lambda = data.InitialLoadFactor;          
%     Uhis = zeros(3*nn,MaxIcr);        
%     FreeDofs = find(sum(Ia,2));
    lmd = 0; icrm = 0; MUL = [U,U];
    Fhis = zeros(MaxIcr,1);
    data_out.t_out=zeros(ne,MaxIcr);        %output member force
    data_out.l_out=zeros(ne,MaxIcr);                % member length
    data_out.q_out=zeros(3*nn,MaxIcr);                % generalized coordinate
%     data_out.Kt_aa_out=cell(1,MaxIcr);       %tangent stiffness of truss
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

        fprintf('icrm = %d, lambda = %6.4f\n',icrm,lmd);

        while err>tol && iter<MaxIter
            iter = iter+1;
%% equilibrium & tangent stiffness matrix

            q=[Ia';Ib']\[qa;qb];
            N=reshape(q(1:3*nn),3,[]);
            sld=q(3*nn+1:end);
            l=sqrt(sum((reshape(q,3,[])*C').^2))'; %bar length
%             l_c=S*l;

            % member force (of truss)
            %         q=E.*A.*(1./l0-1./l);      %force density
            l0s=l0-C*sld;
            strain=(l-l0s)./l0s;        %strain of member
            [E,sigma]=stress_strain(consti_data,index_b,index_s,strain,material);
            t=sigma.*A;         %member force

           
            % equilibrium matrix (of truss)
            H=N*C';
            Cell_H=mat2cell(H,3,ones(1,size(H,2)));          % transfer matrix H into a cell: Cell_H
            A_2=kron(C',eye(3))*blkdiag(Cell_H{:})*diag(l.^-1);     % equilibrium matrix
            Kn=kron(C'*diag(l.\t)*C,eye(3));       %stiffness matrix of truss

            % unbalanced force
            IFa=E_qa'*(-[K_n*N(:);C'*t]+w);

%             IF=w-(Ki*q+phTpn*M);                   %unbalanced force
%             IFa=Ia'*(w-(Ki*q+phTpn*M));            %unbalanced force
%             Fp_a=Ia'*Fp;                   %see the norm of unbalanced force
            % tangent stiffness matrix (of truss)
            K_T=[Kn+A_2*diag(E.*A./l0)*A_2'-A_2*diag(l.\t)*A_2',A_2*diag(E.*A./l0)*C;...
            C'*diag(E.*A./l0)*A_2',C'*diag(E.*A./l0)*C];
            K_T=0.5*(K_T+K_T');
            K_Taa=E_qa'*K_T*E_qa;

            %modify the stiffness matrix
        [V_mode,D]=eig(K_Taa);                       %刚度矩阵特征根
        d=diag(D);                            %eigen value
        lmd=min(d);                     %刚度矩阵最小特征根
        if lmd>0
            Km=K_Taa+u*eye(size(K_Taa)); %修正的刚度矩阵
        else
            Km=K_Taa+(abs(lmd)+u)*eye(size(K_Taa));
        end

           
            dqa = K_Taa\IFa;
        x=1;
        % line search
        if use_energy==1
            opt=optimset('TolX',1e-5);
            [x,V]=fminbnd(@energy_CTS,0,1e1,opt);
        end
        Xa=Xa+x*dXa;


            err = norm(IFa);
            qa=qa+dqa;
            fprintf(' iter = %d, err = %6.4f\n',iter,err);
            if err > 1e8, disp('Divergence!'); break; end
        end

    data_out.q_out(:,icrm)=q;
    data_out.l_out(:,icrm)=l;
    data_out.t_out(:,icrm)=t;      %member force

    end

























