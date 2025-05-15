function ratio_max=ratio_max(h2,R,p,sgg,gbj,h3);
%h2--上环高度
%R--外环高度
%p--complexity for cable dome
%sgg--三角杆的高
%gbj--小三角杆旁边的角的高度
%h3--下环高度
ratio_max=((h2/R)^2+(cos(pi/p)-sgg*sin(gbj)/R*cos(pi/p))^2)/(2*(h2*h3/(R^2)+(cos(pi/p))^2-sgg*sin(gbj)*cos(pi/p)/(R*cos(pi/6))));
end