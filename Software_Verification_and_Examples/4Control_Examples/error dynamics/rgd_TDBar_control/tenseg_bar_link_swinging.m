% function [] = tenseg_link_rotation()
% af2_x1 = [0.0042 0.0198 0.0545 0.1153 0.2057 0.3277 0.4897 0.7021 0.9357]';
% af2x = af2_x1;
% x_crd = [0;af2_x1];
% x_crd = [x_crd;1];
clc;clear all; close all; 
load('member_info.mat');
tenseg_plot(N,CBT,CST); grid on;
%%
qq = length(N(1,:));
N_mid_foil = N(:,1:(qq-1)/3+1);
N_up_foil = N(:,(qq-1)/3+2:2*(qq-1)/3+1);
N_down_foil = N(:,2*(qq-1)/3+2:end);
x_crd =N_mid_foil(1,:); y_crd = N_mid_foil(2,:);
% x_crd = 0:0.1:1;
x_coord=x_crd; y_coord=y_crd;
%% EOM
lq=zeros(1,length(x_coord)-length(x_coord));
for i = 1: length(x_coord)-1
    lq(i) = x_coord(i+1)-x_coord(i);
    %     lq(i) = sqrt((x_coord(i+1)-x_coord(i))^2+(y_coord(i+1)-y_coord(i))^2);
end 
%%
A_p = -.3; % Maximum amplitude
f = 1; % Flapping frequency
EOM= strcat('y=',num2str(A_p),'*sin(2*pi*',num2str(f),'*t)');
tf = 0.3; dt =0.01;
t = 0:dt:tf; % Simulation time
l_r = x_coord(end)-x_coord(1); % tensegrity length
A_m = zeros(1,length(x_coord)-1);
A_m(1)=0;
p1 = [0 0 -1]; p2 = [0 0 1]; % rotation axis
%%
for j = 2:length(x_coord)-1
    lq_sm=0;
    for k=1:j-1
        lq_sm=lq_sm+lq(k);
    end
    A_m(j)=A_p*(lq_sm/l_r)^2;
end
%%
Psi=zeros(1,length(x_coord)-1);
for i = 1:length(x_coord)-2
    Psi(i)=asin((A_m(i+1)-A_m(i))/lq(i));
end
Psi(length(x_coord)-1)= asin((A_p-A_m(end))/lq(end));

y_t=zeros(length(t),length(x_coord)-1);
Psi_t=zeros(length(t),length(x_coord)-1);
for i = 1:length(x_coord)-1
    y_t(:,i)=A_m(i)*sin(2*f*pi.*t);
    Psi_t(:,i)=Psi(i)*sin(2*f*pi.*t);
end

% rotation center
Jx_t=zeros(length(t),length(x_coord));
Jy_t=zeros(length(t),length(x_coord));

% Jx_t(:,1)=x_coord(1);
% Jy_t(:,1)=y_coord(1);


for j=2:length(x_coord)
    for i=1:length(t)
        Jx_t(i,:)=x_coord; Jy_t(i,:)=y_coord;

        Jx_t(i,j) = Jx_t(i,j-1) + lq(j-1)*cos(Psi_t(i,j-1));
        Jy_t(i,j) = Jy_t(i,j-1) + lq(j-1)*sin(Psi_t(i,j-1));
    end
end
C_s_in_add = [];
for i = 1:length(x_coord)-2
    C_s_in_add = [C_s_in_add;length(x_coord)+i length(x_coord)+1+i];
    C_s_in_add = [C_s_in_add;2*(length(x_coord))+i-1 2*(length(x_coord))+i];
    C_s_in_add = [C_s_in_add;i length(x_coord)+1+i];
    C_s_in_add = [C_s_in_add;i 2*(length(x_coord))+i];
end
for i = 1:length(x_coord)-1
    C_s_in_add = [C_s_in_add;i+1 length(x_coord)+i];
    C_s_in_add = [C_s_in_add;i+1 2*(length(x_coord))+i-1];
end
Cb_in = [];
for i = 1: length(Jy_t(1,:))-1
    Cb_in = [Cb_in ; i i+1];
end
Cb_in_add = [];
for k = 1: length(x_coord)-1
    Cb_in_add = [Cb_in_add ; k length(x_coord)+k];
    Cb_in_add = [Cb_in_add ; k 2*length(x_coord)+k-1];
end
Cb_in = [Cb_in; Cb_in_add];

figure(1)
% set(gcf,'Position',get(0,'ScreenSize'));
for i = 1:length(Jx_t(:,1))
    N = [Jx_t(i,1:end);Jy_t(i,1:end);zeros(size(Jx_t(i,1:end)))];
    N_up =  [N(1,1:end-1);N_up_foil(2,:);0*ones(1,length(x_coord)-1)];
    N_down = [N(1,1:end-1);N_down_foil(2,:);0*ones(1,length(x_coord)-1)];
    % rotate
    angle = Psi_t(i,:);
    for j = 1: length(x_coord)-2
        N_up(:,j+1) = trans_space_to_rota((N_up(:,j+1)-N(:,j+1))',p1,p2,angle(j));
        N_up(:,j+1) = N_up(:,j+1) + N(:,j+1);
        N_down(:,j+1) = trans_space_to_rota((N_down(:,j+1)-N(:,j+1))',p1,p2,angle(j));
        N_down(:,j+1) = N_down(:,j+1) + N(:,j+1);
    end
    N_new = [N N_up N_down];
    C_s=tenseg_ind2C(C_s_in_add,N_new);
    C_b = tenseg_ind2C(Cb_in,N_new);
    tenseg_plot(N_new,C_b,C_s,1); grid on;
%     xlim([0 1]);ylim([-.4 .4]); hold on;
    xlabel('X (m)','Interpreter','latex'); ylabel('Y (m)','Interpreter','latex');
    title([EOM,newline,'Real Time:',' ', num2str(i*dt),'s'],'Interpreter','latex');
    tenseg_savegif('bar_link_swinging_test');
    hold off;
    pause(0.1);
end
% beam_swinging_kinematics
% bar_link_swinging_test
