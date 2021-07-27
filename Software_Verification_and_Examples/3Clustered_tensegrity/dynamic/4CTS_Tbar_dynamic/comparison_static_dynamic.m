%%%%% this file plot the dynamic and static comparsion of T bar 
clc; clear all; close all;
%% static result
load cable_net_CTS_static.mat
t_t_sta=t_t;
n_t_sta=n_t;

% plot member force 
tenseg_plot_result(1:substep,t_t_sta([1,2,3,5,6],:),{'1','2','3-4','5','6'},{'Load step','Force (N)'},'plot_member_force.png',saveimg);
grid on;
% Plot nodal coordinate curve X Y
tenseg_plot_result(1:substep,n_t_sta([3*2-1],:),{'2Y'},{'Substep','Coordinate (m)'},'plot_coordinate.png',saveimg);
grid on;

%% dynamic rusult 1s
load cable_net_CTS_dynamic1.mat
t_t_dyn1=t_t;
n_t_dyn1=n_t;
out_tspan_1=out_tspan;

% plot member force 
tenseg_plot_result(out_tspan_1,t_t_dyn1([1:5],:),{'1','2','3','4','5'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);

% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan_1,n_t_dyn1([3*2-1],:),{'2Y'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
hold on
plot(linspace(0,1,substep),n_t_sta([3*2-1],:),'--','linewidth',2);
legend('Dynamic solution','Quasi-static solution');
