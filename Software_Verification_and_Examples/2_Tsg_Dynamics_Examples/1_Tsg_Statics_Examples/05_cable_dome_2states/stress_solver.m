function [sigma_gp_te_error] = stress_solver(A_gp_te)
%this function output stress of top element given the cross sectional area
%of top elements
global A_gp_be I_be I_te data sigmas Gp c_s t l A_gp
%% renew cross sectional area
A_gp=[I_be,I_te]*[A_gp_be;A_gp_te];
A=Gp*A_gp;
%% renew rest length
E=data.E;
l0=E.*A.*l./(t+E.*A);
%% equilibrium calculation
% input data
data.A=A;
data.l0_t=l0;
% nonlinear analysis
data_out=static_solver(data);        %solve equilibrium using mNewton method
% data_out=static_solver2(data);        %solve equilibrium using mNewton method
% data_out{i}=equilibrium_solver_pinv(data);        %solve equilibrium using mNewton method

t_t=data_out.t_out;          % member force in every step
t_gp2=pinv(Gp)*t_t;         % member force in group
sigma_gp=t_gp2./A_gp;       % stress level in group
sigma_gp_be=I_be'*sigma_gp;   %stress level of bottom members

sigma_gp_te=I_te'*sigma_gp;   %stress level of top members
sigma_gp_te_error=sigma_gp_te-sigmas*c_s;
norm(sigma_gp_te_error)
end

