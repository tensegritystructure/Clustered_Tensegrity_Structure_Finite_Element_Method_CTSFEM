function X = skew(x)
% this function is to calculate 3*3 skew symmetric matrix from a 3*1 vector
% or to calculate 3nn_q*3 matrix from 3nn_1*1 vector
x_1=reshape(x,3,[]);

X=zeros(numel(x),3);
for i=1:numel(x)/3;
    X(3*i-2:3*i,:)=[0 -x_1(3,i) x_1(2,i) ; x_1(3,i) 0 -x_1(1,i) ; -x_1(2,i) x_1(1,i) 0 ];
end



% X=[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
end

