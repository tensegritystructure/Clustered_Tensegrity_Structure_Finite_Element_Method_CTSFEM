function V=energy_RDT(x)
% calculate the total energy, x is the cofficient in line search
global   A E l0 E_qa E_qb C qa qb dqa  w nn

q=E_qa*(qa+x*dqa)+E_qb*qb;
 N=reshape(q(1:3*nn),3,[]);
l=sqrt(sum((N*C').^2))'; %bar length 

V=0.5*(l-l0)'*diag(E.*A./l0)*(l-l0)-w'*q;  %结构总势能
