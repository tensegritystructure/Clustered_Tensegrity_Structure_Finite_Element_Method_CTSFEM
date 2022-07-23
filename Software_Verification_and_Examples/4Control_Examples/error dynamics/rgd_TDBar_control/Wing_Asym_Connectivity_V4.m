% Wing Topology 1
% Given: complexity q, wing dimensions
% Solve: Node matrix, connectivity matrices

function [N,CBT,CST] = Wing_Asym_Connectivity_V4(iaf,q,NC,NU,NL)

    nn = 2+3*q; %Number of nodes
    nb = 3*q+1; %Number of bars
    ns = 4*q; %Number of strings
N=zeros(3,nn);
Cb_in=zeros(nb,2);
Cs_in=zeros(ns,2);

% Coordinates calculation

% x_coord = x_coord;

% x_coord=zeros(q+2,1);
% for j=1:q+2
% x_coord(j)=1/(q+1)*(j-1);
% end

% y_coord=zeros(q,1);
% for j=1:q
%     YY=YcoordinateArrayUL(iaf,x_coord(j+1));
% y_coord(j)=YY(2,1);
% end

%% Node matrix
for i=1:q+2
    N(:,i)=NC(:,i);
end

for i=q+3:2*q+2
    N(:,i)=NU(:,i-(q+1));
end

for i=2*q+3:3*q+2
    N(:,i)=NL(:,i-(2*q+1));
end

%% Cb_in

for i=1:q+1
    Cb_in(i,:)=[i,i+1];
end
for i=q+2:2*q+1
     Cb_in(i,:)=[i-q,i+1];
end
for i=2*q+2:3*q+1
     Cb_in(i,:)=[i+1,i-2*q];
end

%% Cs_in
% Upper
Cs_in(1,:)=[1,q+3]; Cs_in(q+1,:)=[2*q+2,q+2]; % 上侧较为特殊

for i=2:q
    Cs_in(i,:)=[i+1+q,i+2+q];
end

%Upper_Inner
for i=q+2:2*q
     Cs_in(i,:)=[i-q,i+2];
end

%Lower
Cs_in(2*q+1,:)=[1,2*q+3]; Cs_in(3*q+1,:)=[3*q+2,q+2];% 下侧较为特殊

for i=2*q+2:3*q
    Cs_in(i,:)=[i+1,i+2];
end

% Lower_Inner
for i=3*q+2:4*q
     Cs_in(i,:)=[i-3*q,i-q+2];
end

CBT=tenseg_ind2C(Cb_in,N);
CST=tenseg_ind2C(Cs_in,N);
%Add a 'rigid'function to update bars and strings
% if RGD==0
% CBT=tenseg_ind2C(Cb_in,N);
% CST=tenseg_ind2C(Cs_in,N);
% else
% [Cb_in_new,Cs_in_new]=Rigid_Convert(RGD,Cb_in,Cs_in,q);
% CBT=tenseg_ind2C(Cb_in_new,N);
% CST=tenseg_ind2C(Cs_in_new,N);
% end


end