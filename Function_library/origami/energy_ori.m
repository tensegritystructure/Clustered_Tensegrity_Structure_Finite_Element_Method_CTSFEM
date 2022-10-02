function V=energy_ori(x)
% calculate the total energy, x is the cofficient in line search
global   A E l0 Ia Ib C Xa Xb dXa  w S theta_0 k_h E_n_total node_in_hinge

X=Ia*(Xa+x*dXa)+Ib*Xb;
N=reshape(X,3,[]);      % nodal coordinate matrix
l=sqrt(sum((reshape(X,3,[])*C').^2))'; %bar length 
l_c=S*l;
V_truss=0.5*(l_c-l0)'*diag(E.*A./l0)*(l_c-l0);     % energy in truss
%% hinge
[~,~,theta]=jacobian_ori(node_in_hinge,N,E_n_total);       % angle



V_hinge=0.5*(theta-theta_0)'*diag(k_h)*(theta-theta_0); % energy in hinge

V=V_truss+V_hinge-w'*X;% total energy
