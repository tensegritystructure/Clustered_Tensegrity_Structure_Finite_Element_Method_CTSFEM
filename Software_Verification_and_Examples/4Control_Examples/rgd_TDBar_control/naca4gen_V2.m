function [iaf] = naca4gen_V2(max_camber,maxcamber_position,thickness,Cosine_Linear,TE)
NACA=[num2str(max_camber), num2str(maxcamber_position), num2str(thickness)];
iaf.designation=NACA;
iaf.HalfCosineSpacing=Cosine_Linear;
iaf.wantFile = 0;
iaf.datFilePath='./'; % Current folder
iaf.is_finiteTE=TE;
% af = naca4gen(iaf);


% iaf1.designation='0012';
% iaf1.n=length(af2_x1);
% iaf1.HalfCosineSpacing=1;
% iaf1.wantFile=1;
% iaf1.datFilePath='./'; % Current folder
% iaf1.is_finiteTE=0;
% af1 = naca4gen(iaf1);
% plot(af1.x,af1.z,'-r*')
% axis equal
% hold on
% A1 = polyarea(af1.x,af1.z);
% 
% A=eye(length(af2_x1));
% B=[1;1];
% af2_x2=kron(A,B);
% % af2_x2 = [1 0 0 0;
% %           1 0 0 0;
% %           0 1 0 0;
% %           0 1 0 0;
% %           0 0 1 0;
% %           0 0 1 0;
% %           0 0 0 1;
% %           0 0 0 1];
% af2_x = af2_x2*af2_x1;
% af2_x = [0 ; af2_x];
% af2_x = [af2_x ; 1];
% af2_z = YcoordinateArrayUL(iaf,af2_x1);
% af2_z = [0  ; af2_z];
% af2_z = [af2_z ; 0];
% 
% af2_xx=zeros(size(af2_x));
% af2_xx(1)=0;
% af2_xx(2:end/2) = af2_x(2:2:end-1);
% af2_xx(end/2+1)=1;
% af2_xx(end/2+2:end) = flipud(af2_x(3:2:end-1));
% af2_xx(end+1)=0;
% af2_zz=zeros(size(af2_z));
% af2_zz(1)=0;
% af2_zz(2:end/2) = af2_z(2:2:end-1);
% af2_zz(end/2+1)=0;
% af2_zz(end/2+2:end) = flipud(af2_z(3:2:end-1));
% af2_zz(end+1)=0;
% 
% 
% % af2 = sturct('x',af2_xx,'z',af2_zz);
% af2.x=af2_xx;
% af2.z=af2_zz;
% % Ae = polyarea(af2_xx(1:end),af2_zz(1:end));
% % plot(af2_xx,af2_zz,'-bo');
% % axis equal
% % legend('Continuous Shape','Cosine Spacing','Error Bound Spacing')
% % legend('Continuous Shape','Error Bound Spacing')
% 
% % xlabel('X-Corrdinate along the cord')
% % ylabel('Y-Coordinate')
% % title('Cosine Spacing v.s. Error Bound Spacing')
% 
% % Ae/A1
% % 
% % e1 = A1/A0
% % e2 = Ae/A0

end