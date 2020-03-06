%------------------------------------------------------------
%     This code replicates Hugget 1993 Environment
%----------	  Written by Suat AKBULUT   ---------------------
%------------   Istanbul, April 2015   ----------------------
%------------------------------------------------------------

%------------------------------------------------------------
%---------------------House Keeping--------------------------
%------------------------------------------------------------
clear;
close all;
clc;

%------------------------------------------------------------
%---------------------Parameters-----------------------------
%------------------------------------------------------------
beta    = 0.96;                                                                % discount rate
alpha   = 1.5;                                                                 % curvature of the utility function
amin    = -2;                                                                  % amin represents the borrowing constraint 
ainc    = 0.01;                                                                %
amax    = 5;                                                                   %
a       = (amin:ainc:amax);                                                    % a grid is the possible values for net asset holding position of a household
na      = size(a,2);
y       = [1 0.1];                                                             % income generated when employed and unemployed, respectively. y(1)=y(e), y(2)=y(u)
sprob   = [0.80 0.20; 0.50 0.50];                                              % s11=prob{s'=e|s=e}, s12=prob{s'=u|s=e}, s21=prob{s'=e|s=u}, s22=prob{s'=u|s=u}
ns      = size(sprob,2);
tol     = 0.00000001;                                                          % tolerance level for VIF iteration
A0      = 1;                                                                   % initial value for aprimenet
suat    = 1;
q       = 0.98                                                                 % Initial guess for q
upper   = 1;                                                                   % Initial upper level for interval including q
lower   = beta;                                                                % Initial lower level for interval including q

aprimenet=A0;
while abs(aprimenet) > 0.0001
    
%------------------------------------------------------------
%--------------------VIF Iteration---------------------------
%------------------------------------------------------------
 
c=zeros(na,ns,na);                                                           % to construct the consumption matrix
for anext = 1:na;                                                            % 1st dim for a'
    for  snow = 1:ns;                                                        % 2nd dim for s
        for anow = 1:na;                                                     % 3rd dim for a                                                                       
       c(anext,snow,anow)= a(anow) + y(snow) - q*a(anext);                   % snow: 1 for employed 2 fo unemployed
       end
    end
end

c(c<=0)=NaN;                                                                 % to mark negative values                        
utility =(c.^(1-alpha))/(1-alpha);
utility(isnan(utility))=-Inf;                                                % to control the nonnegativity of the consumption
v=zeros(ns,na);                                                              % Initial value for Value Function
apolicy=zeros(ns,na);                                                        % To store the location of a', this will determine my policy fnc a'=g(a)
viter=zeros(na,ns,na);                                                       % All possible values for V', during the process this will maximized over a'

diff=1;
iter=0;

while diff > tol
    for anext = 1:na;
        for snow = 1:ns;
             for anow = 1:na;
            
viter(anext, snow, anow) = utility(anext,snow,anow)+beta*sprob(snow,1)*v(1,anext) + beta*sprob(snow,2)*v(2,anext);
             end
        end
    end
    
    for snow = 1:ns;
        for anow = 1:na;
        
            [vnext(snow,anow), apolicy(snow,anow)] = max(viter(:, snow, anow),[],1);
        end
    end
   
    diff = max(max(abs(vnext-v))); 
    bl=(beta/(1-beta))*(min(min(vnext-v)));                                  % McQueen-Porteus Algorithm
    bh=(beta/(1-beta))*(max(max(vnext-v)));     
    v=vnext+(bh+bl)/2;     
    iter=iter+1;    
end
%------------------------------------------------------------
%----------------Stationary Mu Dist--------------------------
%------------------------------------------------------------

    munow=ones(ns,na)/(ns*na);
    distdiff=1;
   while distdiff > tol
    munext=zeros(ns,na);
    for i=1:ns
        for j=1:na
            for k=1:ns
            munext(k,apolicy(i,j)) =  munext(k,apolicy(i,j)) + sprob(i,k)*munow(i,j);
            end
        end
    end
    distdiff=abs(max(max(munext-munow)));
    munow=munext;
   end
   
%------------------------------------------------------------
%-----------Net Asset holding in the Economy-----------------
%------------------------------------------------------------
   
aprimenet=0;                                                                 %to calculate the net asset position of the economy
for snow=1:ns
    for anext=1:na
        aprimenet= aprimenet + a(apolicy(snow,anext))*munext(snow,anext);    % simply the summation of asset positions of each individiaul
    end
end

%------------------------------------------------------------
%----------------Updating q----------------------------------
%------------------------------------------------------------

    if aprimenet<0          %Binary Search Model- Keep contracting the interval that q shoul fall in
        upper=q;
        q=(q+lower)/2
    else
        lower=q;
        q=(q+upper)/2
    end
    
suat=suat+1;
end

disp('PARAMETER VALUES');
disp('');
disp('    beta      alpha      tol     #iterations  '); 
disp([ beta      alpha      tol     suat]);
disp(''); 
disp('EQUILIBRIUM RESULTS ');
disp('');
disp('  Price       Interest rate        A0    ');
disp([   q  1/q-1  aprimenet]);

%------------------------------------------------------------
%----------------Policy Rule---------------------------------
%------------------------------------------------------------

figure
plot(a,a(apolicy))
line([-2,5],[-2,5], 'color', [1,0,0]);
title('Decision Rules')
xlabel('asset levels')
ylabel('asset choices for the next period')
legend('asset choice of employed','asset choice of unemployed')
legend('boxoff')

%------------------------------------------------------------
%---------------Value Function-------------------------------
%------------------------------------------------------------

figure
plot(a,vnext(1,:),'b');
hold on
plot(a,vnext(2,:),'m');
title('Value Functions')
xlabel('asset levels')
ylabel('Lifetime Max. Utility')
legend('employed','unemployed')
legend('boxoff')

%------------------------------------------------------------
%----------------Distribution--------------------------------
%------------------------------------------------------------

figure
plot(a,munow(1,:),'r');
hold on
plot(a,munow(2,:),'b');
title('Stationary Distribution')
xlabel('asset levels')
ylabel('mass of individuals')
legend('employed','unemployed')
legend('boxoff')
%------------------------------------------------------------
%----------------------Welfare Gain--------------------------
%------------------------------------------------------------
pitilda=sprob^100000;
trans=pitilda(1,:);
trans=trans';
cbar=y*trans;
wf= (cbar^(1-alpha))/(1-beta)*(1-alpha);
wi=sum(sum(vnext.*munow));
gain=wf-wi;
disp('    WF      WI      Welfare Gain  '); 
disp([ wf wi gain]);

%------------------------------------------------------------
%----------------Gini Index and Lorenz Curve-----------------
%------------------------------------------------------------
%Tha last command ([g,l,a]=gini(pop,wealth,true)) should be run after
%closing the other figures. I couldn't figure out why it happens.
incy=[1 0; 0 0.1];
w=incy*munow + munow;
wealth=reshape(w,[1402,1]);
pop=reshape(munow,[1402,1]);
[g,l,a]=gini(pop,wealth,true);
