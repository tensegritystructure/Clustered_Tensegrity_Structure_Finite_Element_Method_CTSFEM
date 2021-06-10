
function F_out=Solve_Xcoords_V4(y,ystart,iaf,Tao)

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

%% Function expression
% Nodes coords on three curves: upper(XXU/YYU), new camber(XXC/YYC), lower(XXL/YYL);
    
    if max(ystart,y(2))<=P
    XXU=subs(Xul,m,ystart);XXC=subs(Xtaol,m,y(1));XXL=subs(Xll,m,y(2));
    YYU=subs(Yul,m,ystart);YYC=subs(Ytaol,m,y(1));YYL=subs(Yll,m,y(2));
    DDY1=subs(Dytaol,m,y(1));
    
    F_out=[DDY1*(YYL- YYU)+(XXL-XXU);
       (YYC- YYU)*(XXL-XXC)-(XXC-XXU)*(YYL-YYC)];
    F_out=double(F_out);

    else if min(ystart,y(2))>=P
    XXU=subs(Xur,m,ystart);XXC=subs(Xtaor,m,y(1));XXL=subs(Xlr,m,y(2));
    YYU=subs(Yur,m,ystart);YYC=subs(Ytaor,m,y(1));YYL=subs(Ylr,m,y(2));
    DDY1=subs(Dytaor,m,y(1));
    
    F_out=[DDY1*(YYL- YYU)+(XXL-XXU);
       (YYC- YYU)*(XXL-XXC)-(XXC-XXU)*(YYL-YYC)];
    F_out=double(F_out);   
    
        else if ystart>=y(2) && y(1)>P
    XXU=subs(Xur,m,ystart);XXC=subs(Xtaor,m,y(1));XXL=subs(Xll,m,y(2));
    YYU=subs(Yur,m,ystart);YYC=subs(Ytaor,m,y(1));YYL=subs(Yll,m,y(2));
    DDY1=subs(Dytaor,m,y(1));
    
    F_out=[DDY1*(YYL- YYU)+(XXL-XXU);
       (YYC- YYU)*(XXL-XXC)-(XXC-XXU)*(YYL-YYC)];
    F_out=double(F_out);         
                      
            else if ystart>=y(2) && y(1)<P
    XXU=subs(Xur,m,ystart);XXC=subs(Xtaor,m,y(1));XXL=subs(Xll,m,y(2));
    YYU=subs(Yur,m,ystart);YYC=subs(Ytaor,m,y(1));YYL=subs(Yll,m,y(2));
    DDY1=subs(Dytaol,m,y(1));
    
    F_out=[DDY1*(YYL- YYU)+(XXL-XXU);
       (YYC- YYU)*(XXL-XXC)-(XXC-XXU)*(YYL-YYC)];
    F_out=double(F_out);                  
                else if ystart<y(2) && y(1)>=P
    XXU=subs(Xul,m,ystart);XXC=subs(Xtaor,m,y(1));XXL=subs(Xlr,m,y(2));
    YYU=subs(Yul,m,ystart);YYC=subs(Ytaor,m,y(1));YYL=subs(Ylr,m,y(2));
    DDY1=subs(Dytaor,m,y(1));
    
    F_out=[DDY1*(YYL- YYU)+(XXL-XXU);
       (YYC- YYU)*(XXL-XXC)-(XXC-XXU)*(YYL-YYC)];
    F_out=double(F_out);                      
                    else 
    XXU=subs(Xul,m,ystart);XXC=subs(Xtaol,m,y(1));XXL=subs(Xlr,m,y(2));
    YYU=subs(Yul,m,ystart);YYC=subs(Ytaol,m,y(1));YYL=subs(Ylr,m,y(2));
    DDY1=subs(Dytaol,m,y(1));
    
    F_out=[DDY1*(YYL- YYU)+(XXL-XXU);
       (YYC- YYU)*(XXL-XXC)-(XXC-XXU)*(YYL-YYC)];
    F_out=double(F_out);              
                    end
                end
            end
        end
    end
end
                


