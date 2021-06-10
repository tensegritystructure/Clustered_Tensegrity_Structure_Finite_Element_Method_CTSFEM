clc;clear all; close all;
load('member_info.mat');
tenseg_plot(N,CBT,CST); grid on; hold on;
S_len = []; S_0 = N*CST'; S_len0 = [];
for ss = 1:length(CST(:,1))
    S_len0 = [S_len0;norm(S_0(:,ss),2)];
end
B_len = []; B_0 = N*CBT'; B_len0 = [];
for bb = 1:length(CBT(:,1))
    B_len0 = [B_len0;norm(B_0(:,bb),2)];
end
%%
foil_str = '2412'; propu = min(N(1,2:end)); propl = 1;
for i=1:length(N(1,:))
    if N(2,i)<0
        propl = min(propl,N(1,i));
    end
end
[xU,xL,zU,zL] = gene_airfoil_front(foil_str,propu,propl);
x_front = [xU;xL]; z_front =[zU;zL];
plot(x_front,z_front,'linewidth',2); fig = gcf;
fill(x_front,z_front,[25/255 25/255 112/255]); hold off;
%%
qq = length(N(1,:)); N_mid_foil = N(:,1:(qq-1)/3+1);
N_up_foil = N(:,(qq-1)/3+2:2*(qq-1)/3+1);
N_down_foil = N(:,2*(qq-1)/3+2:end);
x_crd =N_mid_foil(1,:); y_crd = N_mid_foil(2,:);
p1 = [0 0 -1]; p2 = [0 0 1]; % rotation axis
x_coord=x_crd; y_coord=y_crd; 
alpha = -linspace(0,pi/6,6); % alpha = -ones(1,6)*pi/36;
alpha(1) = [];
figure(2); 
% set(gcf,'Position',get(0,'ScreenSize'));
q = (qq-1)/3; mid_q = N_mid_foil(:,(qq-1)/3+1);
dt = 0.01;  S_error = []; B_error = []; N_target = [];
for  i = 1:1/dt
    angle = dt*alpha;
    for j = 1: length(x_coord)-2
        N(:,[j+1:q+1 j+q+2:2*q+1 2*q+j+2:qq]) = [trans_space_to_rota_V2((N(:,[j+1:q+1 j+q+2:2*q+1 2*q+j+2:qq]))',N_mid_foil(:,j)'+p1,N_mid_foil(:,j)'+p2,angle(j))]';
        N_mid_foil = N(:,1:q+1); N_up_foil = N(:,q+2:2*q+1); N_down_foil = N(:,2*q+2:end);
        N_new = [N_mid_foil N_up_foil N_down_foil];
    end
    N_target = [N_target N_new]; 
    xlim([0 1]);ylim([-.4 .4]); hold on;
    xlabel('X (m)','Interpreter','latex'); ylabel('Y (m)','Interpreter','latex');
    plot(x_front,z_front,'linewidth',2);
    fill(x_front,z_front,[25/255 25/255 112/255]); hold on;
    tenseg_plot(N_new,CBT,CST,2); xlim([0 1]);ylim([-.4 .4]); hold on;
    %     title([EOM,newline,'Real Time:',' ', num2str(i*dt),'s'],'Interpreter','latex');
    tenseg_savegif('foil');
    clf; hold off; %     pause(0.1);
    S_len = []; B_len = []; S = N_new*CST';
    for ss = 1: length(S(1,:))
        S_len = [S_len;norm(S(:,ss),2)];
    end
    S_len_error = abs(S_len - S_len0);
    %     S_len_error = S_len;
    S_error = [S_error S_len_error];
    B = N_new*CBT';
    for bb = 1: length(B(1,:))
        B_len = [B_len;norm(B(:,bb),2)];
    end
    B_len_error = (B_len - B_len0); B_error = [B_error B_len_error];
end
%%
for  i = 1:length(B(1,:))
    plot(1:1/dt,B_error(i,:),'-','linewidth',2);
    set(gca,'fontsize', 12,'linewidth',1.2);
    legend_str{i} = ['Bar ' num2str(i)]; hold on;
    title('Bar Length Error','Interpreter','latex');
end
legend(legend_str);
%%
for  i = 1:length(S(1,:))
    plot(1:1/dt,S_error(i,:),'-','linewidth',2);
    set(gca,'fontsize',12,'linewidth',1.2);
    legend_str{i} = ['String ' num2str(i)]; hold on;
    title('String Length Change','Interpreter','latex');
%     set(gca,'ticklength',1.2*get(gca,'ticklength'));
end
legend(legend_str);

% plot(1:100,B_error(1,:),'--');

%     N_new = [N N_up N_down];
%     C_s=tenseg_ind2C(C_s_in_add,N_new);
%     C_b = tenseg_ind2C(Cb_in,N_new);
%     tenseg_plot(N_new,C_b,C_s,1); grid on;
%     xlim([0 1]);ylim([-.4 .4]); hold on;
%     xlabel('X (m)','Interpreter','latex'); ylabel('Y (m)','Interpreter','latex');
%     title([EOM,newline,'Real Time:',' ', num2str(i*dt),'s'],'Interpreter','latex');
%     tenseg_savegif('bar_link_swinging_test');
%     hold off;
%     pause(0.1);
%
% A_p = -.3; % Maximum amplitude
% f = 1; % Flapping frequency
% EOM= strcat('y=',num2str(A_p),'*sin(2*pi*',num2str(f),'*t)');
% tf = 0.3; dt =0.01;
% t = 0:dt:tf; % Simulation time
% l_r = x_coord(end)-x_coord(1); % tensegrity length
% A_m = zeros(1,length(x_coord)-1);
% A_m(1)=0;
