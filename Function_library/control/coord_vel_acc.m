function [n_ct_t,n_ct_dt,n_ct_ddt]=coord_vel_acc(tspan,n_ct1,n_ct2)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% this function give the control coordinate target in time serise
% Inputs:
% tspan: time step 
% n_ct1: nodal coordinate control target: initial position
% n_ct2: nodal coordinate control target: final position
% Outputs:
% n_ct_t: % nodal coordinate of control target in time history
% n_ct_dt: % velocity of control target in time history
% n_ct_ddt: % acceleration of control target in time history

num=numel(tspan);
M_diff=triu(ones(num))-2*triu(ones(num),1)+triu(ones(num),2);   % matrix to calculate diff of vectors in time history
M_diff(:,1)=0;                                      %(1st term 0)
n_ct_t=n_ct1*ones(1,num)+(n_ct2-n_ct1)*linspace(0,1,num);  % nodal coordinate of control target in time history
n_ct_dt=n_ct_t*M_diff;
n_ct_ddt=n_ct_dt*M_diff;
end

