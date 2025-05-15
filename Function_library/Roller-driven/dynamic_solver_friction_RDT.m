function data_out=dynamic_solver_RDT(data_in)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% This function solve dynamic equations using Runge_Kuta method
%
% Inputs:
%	data_in: data structure describing simulation task
%		[].N: initial node positions
%		[].Ia: transform matrix to get free nodal coordinate
%		[].tspan: time series 
% Outputs:
%	History: data structure containing simulation results
%		[].n_t: %time history of nodal coordinate
%		[].t_t: %time history of members' force
%		[].l_t: %time history of members' length
%% input data

E_qa=data_in.E_qa;
q0=data_in.q;
tspan=data_in.tspan;
q0a_d=data_in.q0a_d;    %initial speed of free coordinates
%% dynamic iteration

% initial value
q0a=E_qa\q0;
% n0a_d=zeros(size(n0a));
Q0a=[q0a;q0a_d];
% Perform simulation
data_out = ode4_RDT(@tenseg_dyn_q_qdot_RDT,tspan,Q0a,data_in);


