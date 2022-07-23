function tenseg_cartoon_gif(Ntrace,N,Cb,Cs,name,target_re_index)

% fig=figure;
% v=VideoWriter(name);
% open(v);
[axis_vec, view_vec] = tenseg_axisview(Ntrace);
%%
foil_str = '2412'; propu = min(N(1,2:end)); propl = 1;
for i=1:length(N(1,:))
    if N(2,i)<0
        propl = min(propl,N(1,i));
    end
end    
[xU,xL,zU,zL] = gene_airfoil_front(foil_str,propu,propl);
x_front = [xU;xL]; z_front =[zU;zL];
% plot(x_front,z_front,'linewidth',2); fig = gcf;
% fill(x_front,z_front,[25/255 25/255 112/255]); hold off;
figure(3); 
% set(gcf,'Position',get(0,'ScreenSize'));
%         xlim([0 1]);ylim([-.4 .4]); hold on;
%     xlabel('X (m)','Interpreter','latex'); ylabel('Y (m)','Interpreter','latex');

for i=1:300:size(Ntrace,3)   
       plot(x_front,z_front,'linewidth',3); hold on;
    fill(x_front,z_front,[25/255 25/255 112/255]); hold on;
%       xlim([0 1]);ylim([-.2 .2]); hold on;
    xlabel('X (m)','Interpreter','latex'); ylabel('Y (m)','Interpreter','latex');
   hold on;
%     tenseg_plot(N_new,CBT,CST,2); xlim([0 1]);ylim([-.4 .4]); hold on;
    %     title([EOM,newline,'Real Time:',' ', num2str(i*dt),'s'],'Interpreter','latex');
%     tenseg_savegif('foil');
%     clf; hold off;  %     pause(0.1);


% tenseg_plot_with_axis(Ntrace(:,:,i),Cb,Cs,[],[],view_vec,axis_vec);
tenseg_plot_target_3D_mov(Ntrace(:,:,i),Cb,Cs,target_re_index,3);   hold on;
tenseg_plot_target_3D(Ntrace(:,:,end),Cb,Cs,target_re_index,3);   hold on;
%  xlim([0 1]);ylim([-.2 .2]); hold on;
    xlabel('X (m)','Interpreter','latex'); ylabel('Y (m)','Interpreter','latex');

tenseg_savegif_forever('example2')
% frame=getframe(gcf);
% writeVideo(v,frame);
% close(gcf)
    clf; hold off; 
end
% close(v)
end
