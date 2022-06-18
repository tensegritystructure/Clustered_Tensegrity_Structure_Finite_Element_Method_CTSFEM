function [t_gp,t]=tenseg_prestress_design2(Gp,l,l_gp,A_2ag,V2,w0a,index_gp,fd)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% This function gives the prestress design of tensegrity, consider boundary
% constraints and group constraints, given deisgned force( vector with
% number of prestress modes entries)  This is a new version, applicable to
% CTS
% 
% Inputs:
%   Gp: group matrix of members
%   l: members' length vector
%   H: members' direction matrix
%	A_1ag: Equilirium matrix with boundary, group constraints, force density as variable.
%   V2: prestress mode
%   w0a: external force in free nodal coordinate
%   index_gp: number of groups with designed force
%   fd: force of designed groups
%
% Outputs:
%	q_gp: force desity vector in group
%	t_gp: force vector in group
%	q: force desity vector 
%	t: force vector 
%% 
I=eye(size(Gp,2));
e_d=I(:,index_gp);        % e_d is the matrix to select group of member with designed force
l_d=e_d'*l_gp;            % length of top center circular strings
% qd=fd./l_d;
z=(e_d'*V2)\(e_d'*(fd-pinv(A_2ag)*w0a));   %self-stress coefficient
%%

t_gp=pinv(A_2ag)*w0a+V2*z;

t=Gp*t_gp;                % force in group
end

