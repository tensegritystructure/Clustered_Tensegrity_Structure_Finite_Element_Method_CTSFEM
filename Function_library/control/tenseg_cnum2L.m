function L=tenseg_cnum2L(N,C_num)
%%%generate the L matrix in control
n=size(N,2); %number of node
num=3*size(C_num{1},1)+size(C_num{2},1)+size(C_num{3},1)+size(C_num{4},1);%number of control variables
ind=[C_num{1}*3-2;C_num{1}*3-1;C_num{1}*3;C_num{2}*3-2;C_num{3}*3-1;C_num{4}*3];   %control index

L=zeros(num,3*n);    %control matrix

for i=1:num
    L(i,ind(i))=1;
end