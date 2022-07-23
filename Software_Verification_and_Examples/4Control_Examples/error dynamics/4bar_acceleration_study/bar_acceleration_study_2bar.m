%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%study of acceleration of a single bar%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=[0 0 0;1 0 0;2 0 0]';
C=[-1 1 0;0 -1 1];
B=N*C';
[ne,nn]=size(C);% ne:No.of element;nn:No.of node

l=sqrt(sum((N*C').^2))'; %length matrix
% Plot the structure to make sure it looks right
tenseg_plot(N,C,[]);
% axis off;
%% Boundary constraints
pinned_X=[]; pinned_Y=[]; pinned_Z=(1:3)';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);

%% mass matrix
mass=[1;1];
M=1/6*kron((abs(C)'*diag(mass)*abs(C)+diag(diag(abs(C)'*diag(mass)*abs(C)))),eye(3));
Maa=Ia'*M*Ia;

F=[0 -1 0 0 0 0 0 0 0]';
acc=Maa\(Ia'*F)