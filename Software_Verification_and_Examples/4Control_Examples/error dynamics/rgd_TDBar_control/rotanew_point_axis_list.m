function P_left = rotanew_point_axis_list(P,A_nor,theta)
P_left = [];
for i=1:length(theta)
    A_nor = A_nor/(sqrt(A_nor*A_nor'));
    Ax = A_nor(1);Px = P(1);
    Ay = A_nor(2);Py = P(2);
    Az = A_nor(3);Pz = P(3);
    % theta = pi*120/180;
    % theta_2 = theta*2;
    axp = [Ay*Pz - Az*Py,Az*Px - Ax*Pz,Ax*Py - Ay*Px];
    P_left_temp= P*cos(theta(i)) + (axp)*sin(theta(i)) + A_nor*(A_nor*P')*(1 - cos(theta(i)));
    P_left = [P_left;P_left_temp];
    % P_right = vpa(P*cos(theta_2) + (axp)*sin(theta_2) + A*(A*P')*(1 - cos(theta_2)),10);
end