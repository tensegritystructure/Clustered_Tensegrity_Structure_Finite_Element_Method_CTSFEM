function M=tenseg_mass_matrix_RDT(mass,C,N,l)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% This function generate mass matrix of tensegrity structure.
% 
% Inputs:
%   mass: use lumped matrix 1-yes,0-no
%   C: connectivity matrix
%   N: 1 for lumped matrix; 0 for consistent matrix
% Outputs:
%	M: mass matrix 
%%
% switch lumped
%     case 1  % lumped matrix
% M=0.5*kron(diag(diag(abs(C)'*diag(mass)*abs(C))),eye(3));     
%     case 0  % consistent matrix
% M=1/6*kron((abs(C)'*diag(mass)*abs(C)+diag(diag(abs(C)'*diag(mass)*abs(C)))),eye(3));
% end


H=N*C';
% Cell_H=mat2cell(H,3,ones(1,size(H,2)));          % transfer matrix H into a cell: Cell_H
nn=size(N,2);

M_n=1/6*kron((abs(C)'*diag(mass)*abs(C)+diag(diag(abs(C)'*diag(mass)*abs(C)))),eye(3));
temp=kron(abs(C)',ones(3,1)).*kron(ones(nn,1),H)*diag(l.\mass)*abs(C);
M_ns=1/6*(temp+kron(eye(nn),ones(3,1)).*temp);
M_s=1/6*(abs(C)'*diag(mass)*abs(C)+diag(diag(abs(C)'*diag(mass)*abs(C))));
M=[M_n,M_ns;M_ns',M_s];

end

