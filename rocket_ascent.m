% Aquest codi ha estat realitzat per:
% -Alexis Leon Delgado
% -David Morante Torra
% -Ferran Rubio Vallhonrat
% -Juan Garrido Moreno

% El codi determina l'evoluci� de velocitat i altura al llarg de l'ascens
% vertical d'un coet

clear;
clc;
close all
format long

% Operador boole� per a la consideraci� de resist�ncia aerodin�mica
considerem_drag=true;

%% C�lculs previs

% % Definici� de constants f�siques
R_asterisc=8.31432; % [N m/(mol K)]

% Definici� de les condicions inicials
z_0=0;
V_0=0;
t_0=0;

% Definici� de dades del problema

% Par�metres constants
p_c=100*101325;
T_c=3500;
gamma=1.25;
MW=16e-3;
m_punt=2000;
S_ref=25;
massa_inicial=540e3;
% C�lcul de la R del gas
R_g=R_asterisc/MW;

% Temps de combusti� 
t_b=120;

% Relacions d'�rees estudiades
vector_relacio_area=[20 40 60 80];

% Definici� de l'interval de la variable independent (variaci� del temps)
delta_t=0.01;

% Definici� del vector de la variable independent
t=[t_0:delta_t:t_b];

% Es determina el nombre d'elements de les matrius de resultats
N=length(t);

% Inicialitzaci� de matrius d'emmagatzemament de dades
V=[V_0 zeros(1,N-1)];
h=[z_0 zeros(1,N-1)];

% C�lculs previs de la tovera

% C�lcul del mass flow parameter
MFP=sqrt(gamma)*(1/((1+0.5*(gamma-1)*1^2)^((gamma+1)/(2*(gamma-1)))));
% C�lcul de l'�rea de la gola
A_t=m_punt*sqrt(R_g*T_c)/(p_c*MFP);

% Inici d'un bucle for que resol l'ascens per a cada relaci� d'�rees
% escollida

for j=1:length(vector_relacio_area)
    
% Selecci� de la relaci� d'�rees corresponent
relacio_area=vector_relacio_area(j);

% C�lcul del Mach de sortida M_e
% Es defineix l'equaci� que es desitja a�llar, tota a�llada una banda
% de manera que a l'altra quedi un zero
equacio_M_e=@(M_e)-relacio_area+(2/(gamma+1))^((gamma+1)/(2*(gamma-1)))*1/M_e*(1+(gamma-1)/2* M_e^2)^((gamma+1)/(2*(gamma-1)));

% A continuaci� es troba el zero de la funci�.
% El segon argument de la seg�ent funci� �s un nombre orientatiu entorn del qual es
% troba el resultat
M_e=fzero(equacio_M_e,4);

% C�lcul del mass flow parameter de sortida
MFP_sortida=sqrt(gamma)*(M_e/((1+0.5*(gamma-1)*M_e^2)^((gamma+1)/(2*(gamma-1)))));
% C�lcul de la pressi� de sortida
p_e=p_c/((1+ (gamma-1)/2*M_e^2 )^(gamma/(gamma-1)));

% Definici� pr�via de variables

% Definici� de l'empenyiment
F=@(C_F) p_c*A_t*C_F;
% Definici� de la resist�ncia aerodin�mica
D=@(V_v,rho_v,CD_v) 0.5*V_v^2*rho_v*S_ref*CD_v;
% Definici� de la massa puntual del llan�ador
massa=@(t_v) (massa_inicial-m_punt*t_v);

% Definici� de les equacions diferencials a partir de les variables anteriors 

% Equaci� diferencial de l'altura vers el temps
F_ht=@(V_v) V_v;
% Equaci� diferencial de la velocitat vers el temps
F_Vt=@(t_v,V_v,rho_v,CD_v,C_F,g_v) ((F(C_F)-D(V_v,rho_v,CD_v)  )/massa(t_v))-g_v;

for i=1:N-1 % Bucle de resoluci� per Runge-Kutta
    % Els coeficients del Runge-Kutta no s'avaluen tots al mateix valor de
    % les variables dependents. Per aquesta ra�, cal obtenir les propietats
    % termof�siques per a cada altura
    
    % Obtenci� de propietats per al primer coeficient, els quals
    % coincideixen amb les propietats a l'estat inicial de l'interval
    [T1,P1,rho1,a1,g1,CD1]=getproperties(h(1,i),V(1,i));
    % Imposem resist�ncia nul�la en cas que calgui
    if considerem_drag==false
    CD1=0;
    end
    
    % An�lisi del despreniment del flux en la tovera
    if p_e<(0.4*P1) % Cas de despreniment
        % C�lcul del Mach de sortida efectiu
        M_e_prima=sqrt(2/(gamma-1) *(  (p_c/(0.4*P1))^((gamma-1)/gamma)  -1    )  );
        % C�lcul del mass flow parameter efectiu
        MFP_prima=sqrt(gamma)*(M_e_prima/((1+0.5*(gamma-1)*M_e_prima^2)^((gamma+1)/(2*(gamma-1)))));
        % C�lcul del coeficient d'empenyiment
        C_F=(2/(gamma+1))^((gamma+1)/(2*(gamma-1)))*(gamma*M_e_prima+1/M_e_prima)/ (sqrt(1+(gamma-1)/2* M_e_prima^2))-P1/p_c *MFP/MFP_prima;       
    
    else % Cas sense despreniment
        % C�lcul del coeficient d'empenyiment
        C_F=(2/(gamma+1))^((gamma+1)/(2*(gamma-1)))*(gamma*M_e+1/M_e)/ (sqrt(1+(gamma-1)/2* M_e^2))-(P1/p_c)*relacio_area;
    end
           
    % C�lcul dels coeficients 1
    k_1ht = F_ht(V(1,i));
    k_1Vt = F_Vt(t(1,i),      V(1,i),     rho1,     CD1,     C_F,     g1);
         
    % Obtenci� de les propietats termof�siques per als coeficients 2
    [T2,P2,rho2,a2,g2,CD2]=getproperties(h(1,i)+delta_t/2*k_1ht,V(1,i)+delta_t/2*k_1Vt);
    if considerem_drag==false
    CD2=0;
    end
    
    % C�lcul dels coeficients 2
    k_2ht = F_ht(V(1,i)+delta_t/2*k_1Vt);
    k_2Vt = F_Vt(t(1,i)+delta_t/2,      V(1,i)+delta_t/2*k_1Vt,     rho2,     CD2,     C_F,     g2);
    
    % Obtenci� de propietats termof�siques per als coeficients 3
    [T3,P3,rho3,a3,g3,CD3]=getproperties(h(1,i)+delta_t/2*k_2ht,V(1,i)+delta_t/2*k_2Vt);
    if considerem_drag==false
    CD3=0;
    end
    
    % C�lcul dels coeficients 3
    k_3ht = F_ht(V(1,i)+delta_t/2*k_2Vt);
    k_3Vt = F_Vt(t(1,i)+delta_t/2,      V(1,i)+delta_t/2*k_2Vt,     rho3,     CD3,     C_F,     g3);
    
    % Obtenci� de propietats termof�siques per als coeficients 4
    [T4,P4,rho4,a4,g4,CD4]=getproperties(h(1,i)+delta_t*k_3ht,V(1,i)+delta_t*k_3Vt);
    if considerem_drag==false
    CD4=0;
    end
    
    % C�lcul dels coeficients 4
    k_4ht = F_ht(V(1,i)+delta_t*k_3Vt);
    k_4Vt = F_Vt(t(1,i)+delta_t,      V(1,i)+delta_t*k_3Vt,     rho4,     CD4,     C_F,     g4);
    
    % C�lcul dels valors posteriors
    h(1,i+1)=       h(1,i)+     delta_t/6*(k_1ht+2*k_2ht+2*k_3ht+k_4ht);
    V(1,i+1)=       V(1,i)+     delta_t/6*(k_1Vt+2*k_2Vt+2*k_3Vt+k_4Vt);
    
    Velocitats(j,i+1)=V(1,i+1);
    
end




fig1=figure(1); % Altura respecte del temps
set(fig1,'Renderer', 'painters', 'Position', [100 100 470 400]);
plot(t,h/1000,'DisplayName',['$A_e/A_t=$ ', num2str(vector_relacio_area(j))])
ylabel('Altitud $h\;\mathrm{(km)}$','Interpreter','latex')
xlabel('Temps $t\;\mathrm{(s)}$','Interpreter','latex')
legend('Location','northwest','Interpreter','latex')
if considerem_drag==false
    title('Evoluci\''o temporal de l''altura (cas sense resist\`encia)','interpreter','latex')
else
    title('Evoluci\''o temporal de l''altura (cas amb resist\`encia)','interpreter','latex') 
end
set(gca,'TickLabelInterpreter','latex')
if(j==1)
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.3;
end
hold on


% % 
fig2=figure(2); % Velocitat respecte del temps
set(fig2,'Renderer', 'painters', 'Position', [700 100 470 400]);
plot(t,V,'DisplayName',['$A_e/A_t=$ ', num2str(vector_relacio_area(j))])
ylabel('Velocitat $V\;\mathrm{(m/s)}$','Interpreter','latex')
xlabel('Temps $t\;\mathrm{(s)}$','Interpreter','latex')
legend('Location','northwest','Interpreter','latex')
if considerem_drag==false
    title('Evoluci\''o temporal de la velocitat (cas sense resist\`encia)','interpreter','latex')
else
    title('Evoluci\''o temporal de la velocitat (cas amb resist\`encia)','interpreter','latex') 
end
set(gca,'TickLabelInterpreter','latex')
if(j==1)
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.3;
end
hold on

end