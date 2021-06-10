function cartoon(Ntrace,Cb,Cs,name,target_re_index)

% fig=figure;
v=VideoWriter(name);
open(v);
[axis_vec, view_vec] = tenseg_axisview(Ntrace);
for i=1:5:size(Ntrace,3)                   
% tenseg_plot_with_axis(Ntrace(:,:,i),Cb,Cs,[],[],view_vec,axis_vec);
tenseg_plot_target_3D(Ntrace(:,:,i),Cb,Cs,target_re_index);

frame=getframe(gcf);
writeVideo(v,frame);
close(gcf)
end

close(v)


end
