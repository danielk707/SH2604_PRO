%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SIMULAZIONE E GRAFICI

format short g;                                                            %formato numerico

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INSERIMENTO PARAMETRI SIMULAZIONE

step_variazione=input('Inserire la variazione percentuale di portata [da 0 a 1]  ');                    %variazione portata totale
tempo_variazione=input('Inserire il tempo in cui avviene il cambiamento [in s]  ');                     %tempo
slope_variazione=step_variazione/tempo_variazione;

variaz=int2str(step_variazione)*100;

step_gamma=step_variazione;
slope_gamma=slope_variazione*gamma;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SIMULZIONE BOL
i=1; 
           

%simulazione frazione
risposta1=3;                                                               %variabili per la scelta della simulazione
boolean_gamma=1;
sim('MODEL');
b_gamma_UP=bol_gamma;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%visualizzazione grafici 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CORE

%potenza
figure(1);
hold on;
plot(b_gamma_UP(:,1),b_gamma_UP(:,2)*10^-6,'b','Linewidth',4);
hold off;
set(gca,'XLim',[50 4000],'FontName','Arial','Fontsize', 15, 'Fontweight','b');
grid on;
xlabel('Time [s]','FontName','Arial','Fontsize', 15, 'Fontweight','b')
ylabel('Power [MW]','FontName','Arial','Fontsize', 15, 'Fontweight','b')

%reattivit?
figure(2);
hold on;
plot(b_gamma_UP(:,1),b_gamma_UP(:,3)*1e5,'b','Linewidth',4);
hold off;
set(gca,'XLim',[50 4000],'FontName','Arial','Fontsize', 15, 'Fontweight','b');
grid on;
xlabel('Time [s]','FontName','Arial','Fontsize', 15, 'Fontweight','b')
ylabel('\rho [pcm]','FontName','Arial','Fontsize', 15, 'Fontweight','b')

%Temperature fuel
figure(3);
hold on;
plot(b_gamma_UP(:,1),b_gamma_UP(:,9),'r','Linewidth',4);
hold off;
set(gca,'XLim',[50 4000],'FontName','Arial','Fontsize', 15, 'Fontweight','b');
grid on;
xlabel('Time [s]','FontName','Arial','Fontsize', 15, 'Fontweight','b');
ylabel('Temperature fuel in [C]','FontName','Arial','Fontsize', 15, 'Fontweight','b');


%Temperature coolant
figure(30);
hold on;
plot(b_gamma_UP(:,1),b_gamma_UP(:,10),'y','Linewidth',4);
hold off;
set(gca,'XLim',[50 4000],'FontName','Arial','Fontsize', 15, 'Fontweight','b');
grid on;
xlabel('Time [s]','FontName','Arial','Fontsize', 15, 'Fontweight','b');
ylabel('Temperatures coolant average [C]','FontName','Arial','Fontsize', 15, 'Fontweight','b');


%Temperature Na outlet from fertile
figure(4);
hold on;
plot(b_gamma_UP(:,1),b_gamma_UP(:,11),'c','Linewidth',4);
hold off;
set(gca,'XLim',[50 4000],'FontName','Arial','Fontsize', 15, 'Fontweight','b');
grid on;
xlabel('Time [s]','FontName','Arial','Fontsize', 15, 'Fontweight','b');
ylabel('Temperature coolant outlet [C]','FontName','Arial','Fontsize', 15, 'Fontweight','b');


%Componenti reattivit? 
figure(15)
hold on;
plot(b_gamma_UP(:,1),b_gamma_UP(:,4)*1e5,'r','Linewidth',3);
plot(b_gamma_UP(:,1),b_gamma_UP(:,5)*1e5,'b','Linewidth',3);
plot(b_gamma_UP(:,1),b_gamma_UP(:,6)*1e5,'m','Linewidth',3);
plot(b_gamma_UP(:,1),b_gamma_UP(:,7)*1e5,'g','Linewidth',3);
plot(b_gamma_UP(:,1),b_gamma_UP(:,8)*1e5,'y','Linewidth',3);
hold off;
set(gca,'XLim',[50 4000], 'FontName','Arial','Fontsize', 15, 'Fontweight','b');
grid on;
xlabel('Time [s]','FontName','Arial','Fontsize', 15, 'Fontweight','b');
ylabel('\rho [pcm]','FontName','Arial','Fontsize', 15, 'Fontweight','b');
legend('Doppler','Axial','Coolant','Radial','External','Location','Best');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

