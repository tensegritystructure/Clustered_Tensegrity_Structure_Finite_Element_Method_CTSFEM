
function [NU,NC,NL,q]=Solve_Allcoords_V4(XD,XE,XF,Tao,iaf)

%**************************** Function **********************************%

% 1. A non-fitting error bound method based on NACA airfoil
% 2. Error Bound Main for Asymmetric Wing
% 3. Use loop on Fsolve for total discrete nodes;

%**************************** Data Flow ********************************%

% delta:    Error bound to be satisfied;
% xstart:   Cannot be zero; it's locally singular for the solver;
% x:        Unkown variables; 2 by 1 matrix: x(1) is the point with largest 
%error, x(2) is the point we want
%***********************************************************************%

%% Obtain curve 
syms m

M=str2num(iaf.designation(1))/100;
P=str2num(iaf.designation(2))/10;  % Denote the position of the max camber, which divides the curve into two
T=str2num(iaf.designation(3:4))/100;

a0= 0.2969;
a1=-0.1260;
a2=-0.3516;
a3= 0.2843;

if iaf.is_finiteTE ==1
    a4=-0.1015; % For finite thick TE
else
    a4=-0.1036;  % For zero thick TE
end

if P==0
    Ycl=0;
    Theta1=0;
else 
    Ycl=M/P^2*(2*P*m-m^2);
    Theta1=atan(2*M/P^2*(P-m));
end

Ycr=M/(1-P)^2*(1-2*P+2*P*m-m^2);
Theta2=atan(2*M/(1-P)^2*(P-m));
Yt=(T/0.2)*(a0*sqrt(m)+a1*m+a2*m^2+a3*m^3+a4*m^4);

% G 
Xul=m-Yt*sin(Theta1);   Xur=m-Yt*sin(Theta2);
Yul=Ycl+Yt*cos(Theta1); Yur=Ycr+Yt*cos(Theta2); 

Xtaol=m-Yt*sin(Theta1)*(Tao-1)/(Tao+1);   Xtaor=m-Yt*sin(Theta2)*(Tao-1)/(Tao+1);
Ytaol=Ycl+Yt*cos(Theta1)*(Tao-1)/(Tao+1); Ytaor=Ycr+Yt*cos(Theta2)*(Tao-1)/(Tao+1);

Xll=m+Yt*sin(Theta1);   Xlr=m+Yt*sin(Theta2);
Yll=Ycl-Yt*cos(Theta1); Ylr=Ycr-Yt*cos(Theta2);

% DG
Dyul=diff(Yul,m)/diff(Xul,m); Dyur=diff(Yur,m)/diff(Xur,m);

Dytaol=diff(Ytaol,m)/diff(Xtaol,m); Dytaor=diff(Ytaor,m)/diff(Xtaor,m);

%%
q=length(XD)-2;
NU=zeros(3,q); NU(3,:)=0;
NC=zeros(3,q); NC(3,:)=0;
NL=zeros(3,q); NL(3,:)=0;
%% Function expression
% Nodes coords on three curves: upper(XXU/YYU), new camber(XXC/YYC), lower(XXL/YYL);
for i=1:length(XD)
    
    if max(XD(i),XF(i))<=P
        NU(1,i)=subs(Xul,m,XD(i));NC(1,i)=subs(Xtaol,m,XE(i));NL(1,i)=subs(Xll,m,XF(i));
        NU(2,i)=subs(Yul,m,XD(i));NC(2,i)=subs(Ytaol,m,XE(i));NL(2,i)=subs(Yll,m,XF(i));
    else if min(XD(i),XF(i))>=P
            NU(1,i)=subs(Xur,m,XD(i));NC(1,i)=subs(Xtaor,m,XE(i));NL(1,i)=subs(Xlr,m,XF(i));
            NU(2,i)=subs(Yur,m,XD(i));NC(2,i)=subs(Ytaor,m,XE(i));NL(2,i)=subs(Ylr,m,XF(i))
        else if XD(i)>=XF(i) && XE(i)>=P
                NU(1,i)=subs(Xur,m,XD(i));NC(1,i)=subs(Xtaor,m,XE(i));NL(1,i)=subs(Xll,m,XF(i));
                NU(2,i)=subs(Yur,m,XD(i));NC(2,i)=subs(Ytaor,m,XE(i));NL(2,i)=subs(Yll,m,XF(i))
            else if XD(i)>=XF(i) && XE(i)<P
                    NU(1,i)=subs(Xur,m,XD(i));NC(1,i)=subs(Xtaol,m,XE(i));NL(1,i)=subs(Xll,m,XF(i));
                    NU(2,i)=subs(Yur,m,XD(i));NC(2,i)=subs(Ytaol,m,XE(i));NL(2,i)=subs(Yll,m,XF(i))
                else if XD(i)<XF(i) && XE(i)>=P
                        NU(1,i)=subs(Xul,m,XD(i));NC(1,i)=subs(Xtaor,m,XE(i));NL(1,i)=subs(Xlr,m,XF(i));
                        NU(2,i)=subs(Yul,m,XD(i));NC(2,i)=subs(Ytaor,m,XE(i));NL(2,i)=subs(Ylr,m,XF(i))
                    else
                        NU(1,i)=subs(Xul,m,XD(i));NC(1,i)=subs(Xtaol,m,XE(i));NL(1,i)=subs(Xlr,m,XF(i));
                        NU(2,i)=subs(Yul,m,XD(i));NC(2,i)=subs(Ytaol,m,XE(i));NL(2,i)=subs(Ylr,m,XF(i))
                        
                    end
                end
            end
        end
    end
end

