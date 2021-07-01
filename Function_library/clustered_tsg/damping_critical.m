function d_cri=damping_critical(rho,E_c,A_c)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function gives the critical damping vector of all members of tensegrity structures
%
% Inputs:
%   rho: material density of members
%   E_c: Young's modulus of members
%   A_c: area of members
%
% Outputs:
%	d_cri:critical damping vector of all members 
d_cri=2/sqrt(3)*diag(sqrt(rho))*diag(A_c)*sqrt(E_c);
end

