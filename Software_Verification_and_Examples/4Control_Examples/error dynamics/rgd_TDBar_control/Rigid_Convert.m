%   Rigid_Convert: change certain strings to bars to satisfy rigid wing 
%leading edge
%   Given: full tensegrity connections, Rigid head percent to whole chord
%(Value: 0<RGD<=complexity+1)
%   Solve: new tensegrity connections

function [Cb_in_new,Cs_in_new] =Rigid_Convert(RGD,Cb_in,Cs_in,q)
ns=4*q; % No. of original strings
nb=3*q+1; % Np. of original bars
ncvt=2*(2*RGD-1); % No. of strings to be converted
CVT=zeros(ncvt,2); % Store all those members to be converted

% if RGD==0
%     Cb_in_new=Cb_in; Cs_in_new=Cs_in;
% else
    if RGD==q+1;
    Cb_in_new=[Cb_in;Cs_in];
    Cs_in_new=[];
    else
    % For a random number, use loop to define every CVT entry and also
    % delete those entries
    
    % Upper, 1s layer strings
    for i=1:RGD 
        CVT(i,:)=Cs_in(i,:);
        Cs_in(i,1)=0;  % Here is a trick: first make the entries zero; later delete them after the loop
    end
    % Upper, 2nd layer strings
    for i=RGD+1:2*RGD-1
        CVT(i,:)=Cs_in(i-RGD+(q+1),:);
        Cs_in(i-RGD+(q+1),1)=0;
    end
    
    % Lower, 1s layer strings
    for i=2*RGD:3*RGD-1
        CVT(i,:)=Cs_in(i-2*RGD+ns/2+1,:);
        Cs_in(i-2*RGD+ns/2+1,1)=0;
    end
    
    % Lower, 2nd layer strings
    for i=3*RGD:4*RGD-2
        CVT(i,:)=Cs_in(i-3*RGD+ns/2+(q+1)+1,:);
        Cs_in(i-3*RGD+ns/2+(q+1)+1,1)=0;
    end
    Cb_in_new=[Cb_in;CVT];
    % Detect those entries maked in the loop, then clear them
    id=Cs_in(:,1)<=0 & Cs_in(:,1)>=0;
    Cs_in(id,:)=[];
    Cs_in_new=Cs_in;  
    
end

end