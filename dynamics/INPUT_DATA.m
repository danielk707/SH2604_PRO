%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                                                    %
%                                                                 DATI CORE ESFR                                                                     %
%                                                                                                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                                                   MODELLO CON
%                                                                  SPHEREPAC FUEL
%                                                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                                                    %
%                                                                        CORE                                                                        %
%                                                                                                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
format short g;                                                            %formato di visualizzazione numerico

configurazione=1;                                                          
                                                                 
FAs_inner  = 18;                                                           %[]  numero di FAs inner
FAs_outer  = 18;                                                           %[]  numero di FAs outer
FAs        = FAs_inner + FAs_outer;                                        %[]  numero di FAs totali
pin_per_FA = 37;                                                           %[]  number of pins per FA
%pinside=10;                                                               %[]  numero di pin su un lato

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GLOBAL DATA INPUT
potenzatot0 = 15e6;                                                        %[W] Potenza nominale TOTALE
power_frac_inner_fissile = FAs_inner/FAs;                                  %[]  Frazione della potenza prodotta dal fissile inner
power_frac_outer_fissile = FAs_outer/FAs;                                  %[]  Frazione della potenza prodotta dal fissile outer
power_frac_fertile = 1-power_frac_inner_fissile-power_frac_inner_fissile;  %[]  Frazione della potenza prodotta dal fertile

potenza_inner_fissile0 = potenzatot0 * power_frac_inner_fissile;           %[W] Potenza nominale del fissile inner
potenza_outer_fissile0 = potenzatot0 * power_frac_outer_fissile;           %[W] Potenza nominale del fissile outer
potenza_fertile0       = potenzatot0 * power_frac_fertile;                 %[W] Potenza nominale del fertile

T_in0  = 380;                                                              %[C] Temperatura TV in ingresso al core nominale
T_out0 = 420;                                                              %[C] Temperatura TV in uscita dal core nominale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DATI GEOMETRICI PINS INPUT
d_ext_clad = 2 * 0.95E-2;                                                  %[m] diametro esterno della barretta
d_int_clad = 2 * 0.885E-2;                                                 %[m] diametro esterno della barretta
d_ext_fuel = 2 * 0.85E-2;                                                  %[m] diametro pellet fissile
d_ext_fert = 2 * 0.85E-2;                                                  %[m] diametro pellet fertile
d_int_fuel = 0.0;                                                          %[m] diametro internal hole fuel
d_int_fert = 0.0;                                                          %[m] diametro internal hole fertile

pitch      = 2.5191E-2;                                                      %[m] pin pitch

h_fertile    = 0.0;                                                         %[m] altezza zona fertile sotto
h_fissile    = 1.2;                                                         %[m] altezza fissile
h_plenum_up  = 50e-3;                                                       %[m] upper plenum
h_plenum_bot = 913e-3;                                                      %[m] lower plenum
h_plug_bot   = 82e-3;                                                       %[m] bot plug
h_plug_top   = 18e-3;                                                       %[m] top plug

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DATI GEOMETRICI PINS CALCOLATI
A_pin = (d_ext_clad/2)^2*pi;                                               %[m2] area occupata dalla barretta
P_pin = d_ext_clad*pi;                                                     %[m] perimetro esterno pin

spessore_clad = (d_ext_clad-d_int_clad)/2;                                 %[m] spessore della guaina
spessore_gap  = (d_int_clad-d_ext_fuel)/2;                                 %[m] spessore gap He
spessore_fuel = (d_ext_fuel-d_int_fuel)/2;                                 %[m] spessore pastiglia MOX
spessore_fert = (d_ext_fert-d_int_fert)/2;

pin      = pin_per_FA*FAs;                                                 %[]  numero pin totali
% pin_fert = pin;
pin_fin  = pin_per_FA*FAs_inner;
pin_fout = pin_per_FA*FAs_outer;

h_tot_pin = h_fissile + h_fertile + h_plenum_up + ...
            h_plenum_bot + h_plug_top + h_plug_bot;                        %[m] altezza totale pin
h_attiva  = h_fissile + h_fertile;                                         %[m] altezza attiva

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DATI SPACERS INPUT
d_spacer     = 1.0e-3;                                                     %[m] diametro di wire wrap
pitch_spacer = 225e-3;                                                     %[m] pitch spacer                                                       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DATI SPACERS CALCOLATI
A_spacer    = (d_spacer/2)^2*pi;                                           %[m^2] area occupata da spacer
P_spacer    = d_spacer*pi;                                                 %[m] perimetro di uno spacer
N_spacer    = h_tot_pin/pitch_spacer;
R_clad_wire = d_ext_clad/2+d_spacer/2;
L_spacer    = N_spacer*sqrt(h_tot_pin^+(2*pi()*R_clad_wire)^2);
V_spacer    = L_spacer*A_spacer;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DATI HEX-CAN INPUT
%t_hex         = 4.5e-3/2;                                                  %[m] spessore hex-can
t_hex = 1.0E-3;
flat2flat_out = 2*8.2395E-2;                                                  %[m] outer wrapper flat-to-flat 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DATI HEX-CAN CALCOLATI
flat2flat_in = flat2flat_out-2*t_hex;                                      %[m] inner wrapper flat-to-flat

l_in_hex     = flat2flat_in/sqrt(3);                                       %[m] lato interno hex-can
l_out_hex    = flat2flat_out/sqrt(3);                                      %[m] lato esterno hex-can

P_in_hex     = 6*l_in_hex;                                                 %[m] perimetro interno hex-can
P_out_hex    = 6*l_out_hex;                                                %[m] perimetro esterno hex-can

A_in_hex     = 3*sqrt(3)*(l_in_hex)^2/2;                                   %[m2] area interna hex-can
A_out_hex    = 3*sqrt(3)*(l_out_hex)^2/2;                                  %[m2] area hex-can se si usa il lato esterno
A_hex        = A_out_hex-A_in_hex;                                         %[m2] area del wrapper

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TERMOVETTORE SODIUM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WHOLE CORE
%
T_coolant0 = (T_out0 + T_in0)/2;                                           %[C] Temperatura TV media nella zona core attivo+fertile

%cp_Na_av = 1000*(1.6582-8.4790*10^-4*(T_coolant0+273.15)+4.4541*10^-7*(T_coolant0+273.15)^2-2992.6*(T_coolant0+273.15)^-2);  %[J/(kg*K)]
cp_Na_av  = 176.2 - 4.923E-2 .* (T_coolant0+273.15) + ...
            1.544E-5 .* (T_coolant0+273.15).^2 - 1.524E6 .* (T_coolant0+273.15).^(-2);
%rho_Na_av=219+275.32*(1-(T_coolant0+273.15)/2503.7)+511.58*(1-(T_coolant0+273.15)/2503.7)^0.5;                             %[kg/m^3] 
rho_Na_av = 11441 - 1.2795 .* (T_coolant0+273.15);

area_coolant_tot = A_in_hex*FAs-pin*A_pin-pin*A_spacer;                    %[m2] Area totale coolant
area_coolant     = area_coolant_tot/pin;                                   %[m2] Area passaggio per pin
V_Na_core        = area_coolant_tot*(h_fissile+h_fertile);                 %[m3] Volume coolant nella zona attiva
P_wet_tot        = FAs*P_in_hex+pin*P_pin+pin*P_spacer*0.95;               %[m] Perimetro bagnato totale
P_wet            = P_wet_tot/pin;                                          %[m] Perimetro bagnato un canale
s                = pitch/d_ext_clad;
%d_h              = 4*area_coolant/P_wet;                                   %[m] Diametro idraulico
d_h              = d_ext_clad * (2*sqrt(3)/pi * s^2 - 1);

gamma    = potenzatot0/(cp_Na_av*(T_out0-T_in0));                          %[kg/s] Portata TV core
Portata0 = gamma;
v_Na_av  = gamma/(area_coolant_tot*rho_Na_av);                             %[m/s] velocita' media piombo
M_Na     = gamma*h_attiva/v_Na_av;                                         %[kg] Massa piombo nella zona del CORE fuel 

%mu_Na_av=exp(-6.4406-0.3958*log(T_coolant0+273.15)+556.835/(T_coolant0+273.15));  %[Pa*s] 
mu_Na_av = 4.55E-4 * exp(1069/(T_coolant0+273.15));
%k_Na_av=124.67-0.11381*(T_coolant0+273.15)+5.5226*10^-5*(T_coolant0+273.15)^2-1.1842*10^-8*(T_coolant0+273.15)^3;  %[W/(m*K)] Conducibilita' termica TV
k_Na_av  = 9.2 + 0.011 * (T_coolant0+273.15);

Re_av = rho_Na_av*v_Na_av*d_h/mu_Na_av;                                    %[] Numero di Reynolds medio
Pr_av = mu_Na_av*cp_Na_av/k_Na_av;                                         %[] Numero di Prandtl medio
Pe_av = Re_av*Pr_av;                                                       %[] Numero di Peclet medio

%Nu_Na_av = 7.55*(pitch/d_ext_clad)-20/(pitch/d_ext_clad)^13+0.041/(pitch/d_ext_clad)^2*Pe_av^(0.56+0.19*(pitch/d_ext_clad));     %[] Numero di Nusselt,correlazione di Ushakov semplificata
%Nu_Na_av = 5.0 + 0.025 * Pe_av^0.8;

Nu_Na_av = 4.0 + 0.16 * s^5 + 0.33 * s^3.8 * (Pe_av/100)^0.86;

h_Na_av        = k_Na_av/d_h*Nu_Na_av;                                     %[W/(m^2*K)] Coefficiente di scambio convettivo
h_cs_av        = pi()*d_ext_clad*h_attiva*pin*h_Na_av;                     %[W/k] Coefficiente globale di scambio termico tra guaina e TV
T_clad_out0_av = potenzatot0/h_cs_av+T_coolant0;                           %[C] temperatura esterna della guaina in condizioni stazionarie

cp_Na = cp_Na_av;
%cp_Pb=cp_Pb_av;
V_Na  = V_Na_core;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DATI COMBUSTIBILE
%Processo iterativo per il calcolo delle temperature stazionarie di centro barretta del combustibile e conseguentemente delle propriet? fisiche 
%del combustibile che sono funzione della temperatura dello stesso; le propriet? del combustibile sono calcolate sulla temperatura media della pellet. 
%Quest'ultima e' calcolata, con l'approssimazione di avere k=costante, come 1/2(Tmax + Tsup), da TODREAS.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WHOLE CORE
%
smear       = 0.9775;
rho_teorica = 14.4*1000;                                                  %[kg/m^3] Densita teorica a 20 gradi
Vol_fuel    = (d_ext_fuel^2-d_int_fuel^2)*h_attiva*pi()/4;                %[m^3] volume occupato dal combutibile all'interno 
M_fuel      = smear*Vol_fuel*rho_teorica*pin;                             %[kg] Massa totale combustibile in zona attiva

%U_tot=2.9826*10^6;
%U_tot = 13037 * pi * d_ext_clad * h_attiva; 
%U_tot = 1E5;
%U_tot   = 1.262E5;
lambda_clad = 11;  % [W/(m*K)]
lambda_gap  = 0.6; % [W/(m*K)]
%U_tot = 2*pi*h_attiva*FAs*pin_per_FA/(1/(13037*d_ext_clad/2) + ... 
%U_tot = 2*pi*h_attiva*FAs*pin_per_FA/(1/(9633*d_ext_clad/2) + ... 
U_tot = 2*pi*h_attiva*FAs*pin_per_FA/(1/(h_Na_av*d_ext_clad/2) + ... 
    1/lambda_clad * log(d_ext_clad/d_int_clad) + 1/lambda_gap * log(d_int_clad/d_ext_fuel) + 1/(2*20));
T_fuel0 = potenzatot0/U_tot + T_coolant0;                               %[C] temperatura interna pastiglia in condizioni stazionarie

%U1=193.238;
%U2=162.8647;
%U3=-104.0014;
%U4=29.2056;
%U5=-1.9507;
%U6=2.6441;
%P1=311.7866;
%P2=19.629;
%P3=-0.752;
%P4=0;
%P5=0;
%P6=7.0131;

%temp=(T_fuel0+273.15)/1000;
%cp_UO2=U1+2*U2*temp+3*U3*temp^2+4*U4*temp^3+5*U5*temp^4-U6*temp^-2;
%cp_PuO2=P1+2*P2*temp+3*P3*temp^2+4*P4*temp^3+5*P5*temp^4-P6*temp^-2;
%cp_fuel=(1-0.2)*cp_UO2+0.2*cp_PuO2;                                        %[J/(kg*K)] Capacit? termica fuel
cp_fuel  = 182.63; % J/(kg*K) @900K

tau_fuel = (U_tot/M_fuel/cp_fuel)^-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       
%DATI NEUTRONICA (parametri cinetici e coefficienti di controreazione)

   %parametri cinetici SERPENT
    %beta_i(1) = 69.00E-05;                                                  %beta 1 []
    %beta_i(2) = 22.90E-05;                                                  %beta 2 []
    %beta_i(3) = 59.90E-05;                                                  %beta 3 []
    %beta_i(4) = 124.30E-05;                                                 %beta 4 []
    %beta_i(5) = 55.30E-05;                                                  %beta 5 []
    %beta_i(6) = 69.70E-05;                                                  %beta 6 []
    %beta_i(7) = 0.0;                                                        %beta 7 []
    %beta_i(8) = 0.0;                                                        %beta 8 []
    beta_i(1) = 5.4E-5;                                                  %beta 1 []
    beta_i(2) = 6.24E-4;                                                  %beta 2 []
    beta_i(3) = 2.24E-4;                                                  %beta 3 []
    beta_i(4) = 5.4E-4;                                                 %beta 4 []
    beta_i(5) = 0.0011;                                                  %beta 5 []
    beta_i(6) = 4.56E-4;                                                  %beta 6 []
    beta_i(7) = 4E-4;                                                        %beta 7 []
    beta_i(8) = 2E-4;                                                        %beta 8 []
    %beta_tot=3.98452E-03 ;                                                %beta tot SERPENT []
    
    betatot = 0;
    for j = 1:8
        betatot = betatot + beta_i(j);                                         %beta []
    end

    %lambda_i(1) = 2.49177E-02;                                               %lambda 2 [1/s]
    %lambda_i(2) = 4.25244E-02;                                               %lambda 3 [1/s]
    %lambda_i(3) = 1.33042E-01;                                               %lambda 4 [1/s]
    %lambda_i(4) = 2.92467E-01;                                               %lambda 5 [1/s]
    %lambda_i(5) = 6.66488E-01;                                               %lambda 6 [1/s]
    %lambda_i(6) = 1.81419E-00;                                               %lambda 7 [1/s]
    %lambda_i(7) = 0.0;                                                       %lambda 8 [1/s]
    %lambda_i(8) = 0.0;                                                       %lambda 9 [1/s]  
    lambda_i(1) = 0.0125;                                               %lambda 2 [1/s]
    lambda_i(2) = 0.0283;                                               %lambda 3 [1/s]
    lambda_i(3) = 0.0425;                                               %lambda 4 [1/s]
    lambda_i(4) = 0.133;                                               %lambda 5 [1/s]
    lambda_i(5) = 0.2925;                                               %lambda 6 [1/s]
    lambda_i(6) = 0.66;                                               %lambda 7 [1/s]
    lambda_i(7) = 1.635;                                                       %lambda 8 [1/s]
    lambda_i(8) = 3.55;                                                       %lambda 9 [1/s]  
    
    Linv = 0.44600E-6;                                                       %L invariante

    %coefficienti reattivita' Janne
    alfa_rad     = -0.3e-5;                                                %[1/K] 
    alfa_ax      = -0.3e-5;                                                %[1/K]              
    alfa_CRs     =  1e-5;
    alfa_coolant = -0.3e-5;                                                %[1/K]
    K_Doppler    = -226e-5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VISUALIZZAZIONE A VIDEO DEI RISULTATI

str=fprintf('DATI INIZIALI SIMULAZIONE');
disp(' ')  
disp(' ')  

disp('DATI GLOBALI')
str=fprintf('Potenza totale [MW]                          % 10.4g\n', potenzatot0/1e6);

disp('TEMPERATURE')
str=fprintf('Temperatura coolant media core [C]           % 10.4g\n', T_coolant0);
str=fprintf('Temperatura coolant ingresso core [C]        % 10.4g\n', T_in0);
str=fprintf('Temperatura coolant uscita core [C]          % 10.4g\n',T_out0);
disp(' ')  

disp('PORTATE')
str=fprintf('Portata totale core [kg/s]                   % 10.5g\n', Portata0);

disp(' ')  
disp('TEMPERATURE STAZIONARIO INNER FISSILE')
str=fprintf('Temperatura coolant media inner fissile [C]               % 10.4g\n', T_coolant0);
str=fprintf('Temperatura stazionaria fuel in inner fissile [C]         % 10.4g\n', T_fuel0);
disp(' ')  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SIMULAZIONI

i=1;
risposta1=1;

step_reat     = 0;
slope_reat    = 0;
boolean_reat  = 0;
slope_reat    = 0;

step_gamma    = 0;
slope_gamma   = 0;
boolean_gamma = 0;
slope_gamma   = 0;

step_Tin      = 0;
slope_Tin     = 0;
boolean_Tin   = 0;
slope_Tin     = 0;

Ci_psi        = 1;
Ci_eta        = 1;
