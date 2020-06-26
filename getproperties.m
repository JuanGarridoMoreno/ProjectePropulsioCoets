function [Temperatura,Pressio,Densitat,Velocitat_so,gravetat,CD]=getproperties(Z,V_coet)

% En cas de que l'altitud geomètrica sigui negativa, aquesta es considera
% nul·la per al càlcul de propietats
if(Z<0)
   Z=0; 
end

% Definició de constants físiques
R_asterisc=8.31432; % [N m/(mol K)]
M_0=28.964420*1e-3; % [kg/mol] 
N_A=6.02257e26*1e-3; % [1/mol]
g_0_geo=9.80665; % [m^2/(s^2 m']
g_0=9.80665; % [m/s^2]
Gamma=g_0_geo/g_0; % [m'/m]
r_0=6.356766e6; % [m]

% Vector d'altituds geopotencials
H_b=[0 11 20 32 47 51 71 80]*1e3; % [m']
% Vector de gradients de temperatura
L_Mb=[-6.5 0.0 1.0 2.8 0.0 -2.8 -2.0]*1e-3; % [K/m']
beta_viscosa=1.458*10^(-6); % [kg /(s m K^1/2]
S=110.4; %[K] Sutherland's constant
gamma_aire=1.4;

% Càlcul de pressions base (P_b) i temperatures base a escala molecular (T_b) de cada capa de l'atmosfera
% en base a ISO 2533:1975 (atmosfera ISA)
T_Mb=[288.15 zeros(1,length(L_Mb)-1)]; % [K]
P_b=[101325.0 zeros(1,length(L_Mb)-1)]; % [Pa]

for i=1:(length(L_Mb)-1)
    T_Mb(i+1)=T_Mb(i)+L_Mb(i)*(H_b(i+1)-H_b(i));
    if L_Mb(i)==0
        P_b(i+1)=P_b(i)*exp((-g_0_geo*M_0*(H_b(i+1)-H_b(i)))/(R_asterisc*T_Mb(i)));
    else
        P_b(i+1)=P_b(i)*(T_Mb(i)/(T_Mb(i)+L_Mb(i)*(H_b(i+1)-H_b(i))))^(g_0_geo*M_0/(R_asterisc*L_Mb(i)));
    end
end

% Càlcul de l'altitud geopotencial H segons Z
H=Z*Gamma*r_0/(r_0+Z);

% Càlcul de temperatura i pressió
for b=1:(length(L_Mb))
    if H>=H_b(b) && H<H_b(b+1) % Es distingeix el tram de l'atmosfera
       Temperatura=T_Mb(b)+L_Mb(b)*(H-H_b(b));
       if L_Mb(b)==0 % Si existeix gradient de temperatura
           Pressio=P_b(b)*exp((-g_0_geo*M_0*(H-H_b(b)))/(R_asterisc*T_Mb(b)));
       else % Si no existeix gradient de temperatura
           Pressio=P_b(b)*(T_Mb(b)/(T_Mb(b)+L_Mb(b)*(H-H_b(b))))^(g_0_geo*M_0/(R_asterisc*L_Mb(b)));
       end
       break
    end
end

% Càlcul de la gravetat segons Z
gravetat=g_0*(r_0/(r_0+Z))^2;

% Càlcul de la velocitat del so segons T_M
Velocitat_so=sqrt(gamma_aire*R_asterisc*Temperatura/M_0);

% Càlcul de la densitat
Densitat=(Pressio*M_0)./(Temperatura*R_asterisc); 

% Càlcul del Mach de vol
Mach=V_coet/Velocitat_so;

% Càlcul del coeficient de resistència aerodinàmica
if Mach<0.6
        CD=0.15;
    elseif ((Mach>0.6)&&(Mach<=1.1))
        CD=-4.32*Mach^3+11.016*Mach^2-8.5536*Mach+2.24952;
    elseif ((Mach>1.1)&&(Mach<=1.3))
        CD=-Mach^2+2.2*Mach-0.79;
    else
        CD=0.16769+0.17636/(sqrt(Mach^2-1));
end

end