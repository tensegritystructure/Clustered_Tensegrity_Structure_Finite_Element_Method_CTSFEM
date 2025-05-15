function [gr] = cable_ring_gr_copy(gr_num,p)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
%   用于创建单各环的组数

%p=6;          %complexity for cable dome(Outer ring node)；2的倍数
% gr_num=12;      %同一个环，组数划分（1，2，3，4，6）；节点数需为其的整数倍数
  matrix_1 = [];
  Gr = num2cell(matrix_1)';  %杆矩阵
  
switch gr_num
    case 1
%         gr={[2*p+1:4*p]';[1:2*p]';[4*p+1:5*p]'};  % 上，下，环；一组
        gr_fz1=[zeros(1,p/6),11*p/6*ones(1,p/6)];  %赋值
        gr_xxs1=[kron(ones(1,p/6),[1,2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs2=[kron(ones(1,p/6),[p/6+1,p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs3=[kron(ones(1,p/6),[2*p/6+1,2*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs4=[kron(ones(1,p/6),[3*p/6+1,3*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs5=[kron(ones(1,p/6),[4*p/6+1,4*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs6=[kron(ones(1,p/6),[5*p/6+1,5*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs7=[kron(ones(1,p/6),[6*p/6+1,6*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs8=[kron(ones(1,p/6),[7*p/6+1,7*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs9=[kron(ones(1,p/6),[8*p/6+1,8*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs10=[kron(ones(1,p/6),[9*p/6+1,9*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs11=[kron(ones(1,p/6),[10*p/6+1,10*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs12=[kron(ones(1,p/6),[11*p/6+1,11*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_1={gr_xxs1;gr_xxs2;gr_xxs3;gr_xxs4;gr_xxs5;gr_xxs6;gr_xxs7;gr_xxs8;gr_xxs9;gr_xxs10;gr_xxs11;gr_xxs12;
            [4*p+1:6*p]';[6*p+1:7*p]';[7*p+1:9*p]';[9*p+1:11*p]'};  % 上，下，环；一组 ；索矩阵
      %  gr_2={[2*p+1:4*p]'};%斜上小三角
        gr = [gr_1;Gr];

%         gr={[4*p+1:5*p]';[2*p+1:4*p]';[1:2*p]'};  % 上，下，环；一组
    case 2
        gr_fz=[zeros(1,p/2),p/2*ones(1,p/2)];  %赋值
        gr_fz1=[zeros(1,p/6),11*p/6*ones(1,p/6)];  %赋值
        gr_xxs1=[kron(ones(1,p/6),[1,2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs2=[kron(ones(1,p/6),[p/6+1,p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs3=[kron(ones(1,p/6),[2*p/6+1,2*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs4=[kron(ones(1,p/6),[3*p/6+1,3*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs5=[kron(ones(1,p/6),[4*p/6+1,4*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs6=[kron(ones(1,p/6),[5*p/6+1,5*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs7=[kron(ones(1,p/6),[6*p/6+1,6*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs8=[kron(ones(1,p/6),[7*p/6+1,7*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs9=[kron(ones(1,p/6),[8*p/6+1,8*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs10=[kron(ones(1,p/6),[9*p/6+1,9*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs11=[kron(ones(1,p/6),[10*p/6+1,10*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs12=[kron(ones(1,p/6),[11*p/6+1,11*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        
        gr_ssxs1=[kron(ones(1,p/2),[4*p+1,4*p+2])+2*kron(0:(p-1)/2,[1,1])+gr_fz]';
        gr_ssxs2=[kron(ones(1,p/2),[9*p/2+1,9*p/2+2])+2*kron(0:(p-1)/2,[1,1])+gr_fz]';
        gr_zsxs1=[kron(ones(1,p/2),[7*p+1,7*p+2])+2*kron(0:(p-1)/2,[1,1])+gr_fz]';
        gr_zsxs2=[kron(ones(1,p/2),[15*p/2+1,15*p/2+2])+2*kron(0:(p-1)/2,[1,1])+gr_fz]';
        gr_zxxs1=[kron(ones(1,p/2),[9*p+1,9*p+2])+2*kron(0:(p-1)/2,[1,1])+gr_fz]';
        gr_zxxs2=[kron(ones(1,p/2),[19*p/2+1,19*p/2+2])+2*kron(0:(p-1)/2,[1,1])+gr_fz]';
        gr_hs1= [6*p+1:13*p/2]';
        gr_hs2=[13*p/2+1:7*p]';
        gr_1={gr_hs1;gr_hs2;gr_ssxs1;gr_ssxs2;gr_xxs1;gr_xxs2;gr_xxs3;gr_xxs4;gr_xxs5;gr_xxs6;gr_xxs7;gr_xxs8;gr_xxs9;gr_xxs10;gr_xxs11;gr_xxs12;
            gr_zsxs1;gr_zsxs2;gr_zxxs1;gr_zxxs2}; %两组; 环，上，下
        gr = [gr_1;Gr];
        
    case 3
        gr_fz=[zeros(1,p/3),2*p/3*ones(1,p/3)];  %赋值
        gr_fz1=[zeros(1,p/6),11*p/6*ones(1,p/6)];  %赋值
        gr_xxs1=[kron(ones(1,p/6),[1,2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs2=[kron(ones(1,p/6),[p/6+1,p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs3=[kron(ones(1,p/6),[2*p/6+1,2*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs4=[kron(ones(1,p/6),[3*p/6+1,3*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs5=[kron(ones(1,p/6),[4*p/6+1,4*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs6=[kron(ones(1,p/6),[5*p/6+1,5*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs7=[kron(ones(1,p/6),[6*p/6+1,6*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs8=[kron(ones(1,p/6),[7*p/6+1,7*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs9=[kron(ones(1,p/6),[8*p/6+1,8*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs10=[kron(ones(1,p/6),[9*p/6+1,9*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs11=[kron(ones(1,p/6),[10*p/6+1,10*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs12=[kron(ones(1,p/6),[11*p/6+1,11*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        
        gr_ssxs1=[kron(ones(1,p/3),[4*p+1,4*p+2])+2*kron(0:(p-1)/3,[1,1])+gr_fz]';
        gr_ssxs2=[kron(ones(1,p/3),[13*p/3+1,13*p/3+2])+2*kron(0:(p-1)/3,[1,1])+gr_fz]';
        gr_ssxs3=[kron(ones(1,p/3),[14*p/3+1,14*p/3+2])+2*kron(0:(p-1)/3,[1,1])+gr_fz]';
        gr_zsxs1=[kron(ones(1,p/3),[7*p+1,7*p+2])+2*kron(0:(p-1)/3,[1,1])+gr_fz]';
        gr_zsxs2=[kron(ones(1,p/3),[22*p/3+1,22*p/3+2])+2*kron(0:(p-1)/3,[1,1])+gr_fz]';
        gr_zsxs3=[kron(ones(1,p/3),[23*p/3+1,23*p/3+2])+2*kron(0:(p-1)/3,[1,1])+gr_fz]';
        gr_zxxs1=[kron(ones(1,p/3),[9*p+1,9*p+2])+2*kron(0:(p-1)/3,[1,1])+gr_fz]';
        gr_zxxs2=[kron(ones(1,p/3),[28*p/3+1,28*p/3+2])+2*kron(0:(p-1)/3,[1,1])+gr_fz]';
        gr_zxxs3=[kron(ones(1,p/3),[29*p/3+1,29*p/3+2])+2*kron(0:(p-1)/3,[1,1])+gr_fz]';
        gr_hs1= [6*p+1:19*p/3]';
        gr_hs2=[19*p/3+1:20*p/3]';
        gr_hs3=[20*p/3+1:7*p]';
        gr_1={gr_hs1;gr_hs2;gr_hs3;gr_xxs1;gr_xxs2;gr_xxs3;gr_xxs4;gr_xxs5;gr_xxs6;gr_xxs7;gr_xxs8;gr_xxs9;gr_xxs10;gr_xxs11;gr_xxs12;
            gr_ssxs1;gr_ssxs2;gr_ssxs3;gr_zsxs1;gr_zsxs2;gr_zsxs3;gr_zxxs1;gr_zxxs2;gr_zxxs3}; %三组;环，上, 下，
        gr = [gr_1;Gr];

    case 4
        gr_fz=[zeros(1,p/4),3*p/4*ones(1,p/4)];  %赋值
        gr_xxs1=[kron(ones(1,p/4),[1,2])+2*kron(0:(p-1)/4,[1,1])+gr_fz]';
        gr_xxs2=[kron(ones(1,p/4),[p/4+1,p/4+2])+2*kron(0:(p-1)/4,[1,1])+gr_fz]';
        gr_xxs3=[kron(ones(1,p/4),[2*p/4+1,2*p/4+2])+2*kron(0:(p-1)/4,[1,1])+gr_fz]';
        gr_xxs4=[kron(ones(1,p/4),[3*p/4+1,3*p/4+2])+2*kron(0:(p-1)/4,[1,1])+gr_fz]';
        gr_sxs1=[kron(ones(1,p/4),[2*p+1,2*p+2])+2*kron(0:(p-1)/4,[1,1])+gr_fz]';
        gr_sxs2=[kron(ones(1,p/4),[9*p/4+1,9*p/4+2])+2*kron(0:(p-1)/4,[1,1])+gr_fz]';
        gr_sxs3=[kron(ones(1,p/4),[10*p/4+1,10*p/4+2])+2*kron(0:(p-1)/4,[1,1])+gr_fz]';
        gr_sxs4=[kron(ones(1,p/4),[11*p/4+1,11*p/4+2])+2*kron(0:(p-1)/4,[1,1])+gr_fz]';
        gr_hs1= [4*p+1:17*p/4]';
        gr_hs2=[17*p/4+1:18*p/4]';
        gr_hs3=[18*p/4+1:19*p/4]';
        gr_hs4=[19*p/4+1:5*p]';
        gr_1={gr_xxs1;gr_xxs2;gr_xxs3;gr_xxs4;gr_hs1;gr_hs2;gr_hs3;gr_hs4;gr_sxs1;gr_sxs2;gr_sxs3;gr_sxs4}; %四组; 下，环，上
        gr = [gr_1;Gr];
        
    case 6
        gr_fz=[zeros(1,p/6),5*p/6*ones(1,p/6)];  %赋值
        gr_fz1=[zeros(1,p/6),11*p/6*ones(1,p/6)];  %赋值
        gr_xxs1=[kron(ones(1,p/6),[1,2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs2=[kron(ones(1,p/6),[p/6+1,p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs3=[kron(ones(1,p/6),[2*p/6+1,2*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs4=[kron(ones(1,p/6),[3*p/6+1,3*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs5=[kron(ones(1,p/6),[4*p/6+1,4*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs6=[kron(ones(1,p/6),[5*p/6+1,5*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs7=[kron(ones(1,p/6),[6*p/6+1,6*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs8=[kron(ones(1,p/6),[7*p/6+1,7*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs9=[kron(ones(1,p/6),[8*p/6+1,8*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs10=[kron(ones(1,p/6),[9*p/6+1,9*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs11=[kron(ones(1,p/6),[10*p/6+1,10*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
        gr_xxs12=[kron(ones(1,p/6),[11*p/6+1,11*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz1]';
       
       
        gr_ssxs1=[kron(ones(1,p/6),[4*p+2,4*p+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_ssxs2=[kron(ones(1,p/6),[26*p/6+1,25*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_ssxs3=[kron(ones(1,p/6),[27*p/6+1,26*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_ssxs4=[kron(ones(1,p/6),[28*p/6+1,27*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_ssxs5=[kron(ones(1,p/6),[29*p/6+1,28*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_ssxs6=[kron(ones(1,p/6),[4*p+1,29*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        
        gr_zsxs1=[kron(ones(1,p/6),[7*p+2,7*p+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_zsxs2=[kron(ones(1,p/6),[44*p/6+1,43*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_zsxs3=[kron(ones(1,p/6),[45*p/6+1,44*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_zsxs4=[kron(ones(1,p/6),[46*p/6+1,45*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_zsxs5=[kron(ones(1,p/6),[47*p/6+1,46*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_zsxs6=[kron(ones(1,p/6),[7*p+1,47*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        
        gr_zxxs1=[kron(ones(1,p/6),[9*p+2,9*p+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_zxxs2=[kron(ones(1,p/6),[56*p/6+1,55*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_zxxs3=[kron(ones(1,p/6),[57*p/6+1,56*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_zxxs4=[kron(ones(1,p/6),[58*p/6+1,57*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_zxxs5=[kron(ones(1,p/6),[59*p/6+1,58*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_zxxs6=[kron(ones(1,p/6),[9*p+1,59*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
         
        gr_hs1= [6*p+1:37*p/6]';
        gr_hs2=[37*p/6+1:38*p/6]';
        gr_hs3=[38*p/6+1:39*p/6]';
        gr_hs4=[39*p/6+1:40*p/6]';
        gr_hs5=[40*p/6+1:41*p/6]';
        gr_hs6=[41*p/6+1:7*p]';
        gr_1={gr_xxs1;gr_xxs2;gr_xxs3;gr_xxs4;gr_xxs5;gr_xxs6;gr_xxs7;gr_xxs8;gr_xxs9;gr_xxs10;gr_xxs11;gr_xxs12;
            gr_ssxs1;gr_ssxs2;gr_ssxs3;gr_ssxs4;gr_ssxs5;gr_ssxs6;
            gr_zsxs1;gr_zsxs2;gr_zsxs3;gr_zsxs4;gr_zsxs5;gr_zsxs6;
            gr_zxxs1;gr_zxxs2;gr_zxxs3;gr_zxxs4;gr_zxxs5;gr_zxxs6;
             gr_hs1;gr_hs2;gr_hs3;gr_hs4;gr_hs5;gr_hs6}; %六组; 下，上，环
         gr = [gr_1;Gr];
         
       case 12  
       gr_fz=[zeros(1,p/6),11*p/6*ones(1,p/6)];  %赋值
        gr_xxs1=[kron(ones(1,p/6),[1,2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_xxs2=[kron(ones(1,p/6),[p/6+1,p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_xxs3=[kron(ones(1,p/6),[2*p/6+1,2*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_xxs4=[kron(ones(1,p/6),[3*p/6+1,3*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_xxs5=[kron(ones(1,p/6),[4*p/6+1,4*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_xxs6=[kron(ones(1,p/6),[5*p/6+1,5*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_xxs7=[kron(ones(1,p/6),[6*p/6+1,6*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_xxs8=[kron(ones(1,p/6),[7*p/6+1,7*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_xxs9=[kron(ones(1,p/6),[8*p/6+1,8*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_xxs10=[kron(ones(1,p/6),[9*p/6+1,9*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_xxs11=[kron(ones(1,p/6),[10*p/6+1,10*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
        gr_xxs12=[kron(ones(1,p/6),[11*p/6+1,11*p/6+2])+2*kron(0:(p-1)/6,[1,1])+gr_fz]';
       
        gr_ssxs1= [4*p+1:25*p/6]';
        gr_ssxs2=[25*p/6+1:26*p/6]';
        gr_ssxs3=[26*p/6+1:27*p/6]';
        gr_ssxs4=[27*p/6+1:28*p/6]';
        gr_ssxs5=[28*p/6+1:29*p/6]';
        gr_ssxs6=[29*p/6+1:5*p]';
        gr_ssxs7= [5*p+1:31*p/6]';
        gr_ssxs8=[31*p/6+1:32*p/6]';
        gr_ssxs9=[32*p/6+1:33*p/6]';
        gr_ssxs10=[33*p/6+1:34*p/6]';
        gr_ssxs11=[34*p/6+1:35*p/6]';
        gr_ssxs12=[35*p/6+1:6*p]';
        
        gr_zsxs1= [7*p+1:43*p/6]';
        gr_zsxs2=[43*p/6+1:44*p/6]';
        gr_zsxs3=[44*p/6+1:45*p/6]';
        gr_zsxs4=[45*p/6+1:46*p/6]';
        gr_zsxs5=[46*p/6+1:47*p/6]';
        gr_zsxs6=[47*p/6+1:8*p]';
        gr_zsxs7= [8*p+1:49*p/6]';
        gr_zsxs8=[49*p/6+1:50*p/6]';
        gr_zsxs9=[50*p/6+1:51*p/6]';
        gr_zsxs10=[51*p/6+1:52*p/6]';
        gr_zsxs11=[52*p/6+1:53*p/6]';
        gr_zsxs12=[53*p/6+1:9*p]';
        
        gr_zxxs1= [9*p+1:55*p/6]';
        gr_zxxs2=[55*p/6+1:56*p/6]';
        gr_zxxs3=[56*p/6+1:57*p/6]';
        gr_zxxs4=[57*p/6+1:58*p/6]';
        gr_zxxs5=[58*p/6+1:59*p/6]';
        gr_zxxs6=[59*p/6+1:10*p]';
        gr_zxxs7= [10*p+1:61*p/6]';
        gr_zxxs8=[61*p/6+1:62*p/6]';
        gr_zxxs9=[62*p/6+1:63*p/6]';
        gr_zxxs10=[63*p/6+1:64*p/6]';
        gr_zxxs11=[64*p/6+1:65*p/6]';
        gr_zxxs12=[65*p/6+1:11*p]';
        
        gr_hs1= [6*p+1:37*p/6]';
        gr_hs2=[37*p/6+1:38*p/6]';
        gr_hs3=[38*p/6+1:39*p/6]';
        gr_hs4=[39*p/6+1:40*p/6]';
        gr_hs5=[40*p/6+1:41*p/6]';
        gr_hs6=[41*p/6+1:7*p]';

        gr_1={gr_xxs1;gr_xxs2;gr_xxs3;gr_xxs4;gr_xxs5;gr_xxs6;gr_xxs7;gr_xxs8;gr_xxs9;gr_xxs10;gr_xxs11;gr_xxs12;
              gr_ssxs1;gr_ssxs2;gr_ssxs3;gr_ssxs4;gr_ssxs5;gr_ssxs6;gr_ssxs7;gr_ssxs8;gr_ssxs9;gr_ssxs10;gr_ssxs11;gr_ssxs12;
              gr_zsxs1;gr_zsxs2;gr_zsxs3;gr_zsxs4;gr_zsxs5;gr_zsxs6;gr_zsxs7;gr_zsxs8;gr_zsxs9;gr_zsxs10;gr_zsxs11;gr_zsxs12;
              gr_zxxs1;gr_zxxs2;gr_zxxs3;gr_zxxs4;gr_zxxs5;gr_zxxs6;gr_zxxs7;gr_zxxs8;gr_zxxs9;gr_zxxs10;gr_zxxs11;gr_zxxs12;
              gr_hs1;gr_hs2;gr_hs3;gr_hs4;gr_hs5;gr_hs6};
 
 
        gr = [Gr;gr_1];

         
end
end

