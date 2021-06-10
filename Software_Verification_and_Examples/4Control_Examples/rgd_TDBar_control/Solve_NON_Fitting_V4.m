
function F_out=Solve_NON_Fitting_V4(x,xstart,delta,iaf)
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

Xll=m+Yt*sin(Theta1);   Xlr=m+Yt*sin(Theta2);
Yll=Ycl-Yt*cos(Theta1); Ylr=Ycr-Yt*cos(Theta2);

% DG
Dyul=diff(Yul,m)/diff(Xul,m); Dyur=diff(Yur,m)/diff(Xur,m);

%% Function expression

if x(2)<=P
    XXstart=subs(Xul,m,xstart);XX1=subs(Xul,m,x(1));XX2=subs(Xul,m,x(2));
    YYstart=subs(Yul,m,xstart);YY1=subs(Yul,m,x(1));YY2=subs(Yul,m,x(2));
    DD1=subs(Dyul,m,x(1));
    
    A=(YY2-YYstart)/(XX2-XXstart);
    B=-1;
    C=YY2-A*XX2;
    
    F_out=[abs(A*XX1+B*YY1+C)/sqrt(A*A+B*B)-delta;
        DD1-A];
    F_out=double(F_out);
    
else if xstart>=P
        XXstart=subs(Xur,m,xstart);XX1=subs(Xur,m,x(1));XX2=subs(Xur,m,x(2));
        YYstart=subs(Yur,m,xstart);YY1=subs(Yur,m,x(1));YY2=subs(Yur,m,x(2));
        DD1=subs(Dyur,m,x(1));
        
        A=(YY2-YYstart)/(XX2-XXstart);
        B=-1;
        C=YY2-A*XX2;
        
        F_out=[abs(A*XX1+B*YY1+C)/sqrt(A*A+B*B)-delta;
            DD1-A];
        F_out=double(F_out);
        
    else if x(1)<=P
            XXstart=subs(Xul,m,xstart);XX1=subs(Xul,m,x(1));XX2=subs(Xur,m,x(2));
            YYstart=subs(Yul,m,xstart);YY1=subs(Yul,m,x(1));YY2=subs(Yur,m,x(2));
            DD1=subs(Dyul,m,x(1));
            
            A=(YY2-YYstart)/(XX2-XXstart);
            B=-1;
            C=YY2-A*XX2;
            
            F_out=[abs(A*XX1+B*YY1+C)/sqrt(A*A+B*B)-delta;
                DD1-A];
            F_out=double(F_out);
            
        else
            XXstart=subs(Xul,m,xstart);XX1=subs(Xur,m,x(1));XX2=subs(Xur,m,x(2));
            YYstart=subs(Yul,m,xstart);YY1=subs(Yur,m,x(1));YY2=subs(Yur,m,x(2));
            DD1=subs(Dyur,m,x(1));
            
            A=(YY2-YYstart)/(XX2-XXstart);
            B=-1;
            C=YY2-A*XX2;
            
            F_out=[abs(A*XX1+B*YY1+C)/sqrt(A*A+B*B)-delta;
                DD1-A];
            F_out=double(F_out);
            
        end
    end
end
end
