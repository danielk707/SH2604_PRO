%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SIMULAZIONE E GRAFICI

format short g;                                                            %formato numerico
step_reat=0;    
step_Tin=0; 
step_gamma=0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INSERIMENTO PARAMETRI SIMULAZIONE

step_reat=input('Insert the reactivity [in pcm] (nominal value = 0) ');                  %variazione reattivita
reat=int2str(step_reat);
tempo_reat=input('Insert the time [in s]  ');                                            %tempo
slope_reat=step_reat/tempo_reat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SIMULZIONE BOL
i=1; 
           
%simulazione reattivit?
risposta1=1;                                                               %variabili per la scelta della simulazione
boolean_reat=1;
sim('MODEL');
b_reat_UP=bol_reat;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%visualizzazione grafici 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%step REATTIVITA'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CORE

%potenza
figure(1);
set(1,'name',['Ramp ' reat ' pcm, Power']);
hold on;
plot(b_reat_UP(:,1),b_reat_UP(:,2)*10^-6,'Linewidth',4);
hold off;
set(gca,'XLim',[50 1000],'FontName','Arial','Fontsize', 15, 'Fontweight','b');
grid on;
xlabel('Time [s]','FontName','Arial','Fontsize', 15, 'Fontweight','b')
ylabel('Power [MW]','FontName','Arial','Fontsize', 15, 'Fontweight','b')

%reattivit?
figure(2);
set(2,'name',['Ramp ' reat ' pcm, Reactivity']);
hold on;
plot(b_reat_UP(:,1),b_reat_UP(:,3)*1e5,'Linewidth',4);
hold off;
set(gca,'XLim',[50 1000],'FontName','Arial','Fontsize', 15, 'Fontweight','b');
grid on;
xlabel('Time [s]','FontName','Arial','Fontsize', 15, 'Fontweight','b')
ylabel('\rho [pcm]','FontName','Arial','Fontsize', 15, 'Fontweight','b')

%Temperature fuel
figure(3);
set(3,'name',['Ramp ' reat ' pcm, Temperature fuel']);
hold on;
plot(b_reat_UP(:,1),b_reat_UP(:,9),'Linewidth',4);
hold off;
set(gca,'XLim',[50 1000],'FontName','Arial','Fontsize', 15, 'Fontweight','b');
grid on;
xlabel('Time [s]','FontName','Arial','Fontsize', 15, 'Fontweight','b');
ylabel('Temperature fuel in [C]','FontName','Arial','Fontsize', 15, 'Fontweight','b');


%Temperature coolant
figure(30);
set(30,'name',['Ramp ' reat ' pcm, Temperature coolant']);
hold on;
plot(b_reat_UP(:,1),b_reat_UP(:,10),'Linewidth',4);
hold off;
set(gca,'XLim',[50 1000],'FontName','Arial','Fontsize', 15, 'Fontweight','b');
grid on;
xlabel('Time [s]','FontName','Arial','Fontsize', 15, 'Fontweight','b');
ylabel('Temperatures coolant average [C]','FontName','Arial','Fontsize', 15, 'Fontweight','b');


%Temperature Na outlet from fertile
figure(4);
set(4,'name',['Ramp ' reat ' pcm, T Na outlet']);
hold on;
plot(b_reat_UP(:,1),b_reat_UP(:,11),'Linewidth',4);
hold off;
set(gca,'XLim',[50 1000],'FontName','Arial','Fontsize', 15, 'Fontweight','b');
grid on;
xlabel('Time [s]','FontName','Arial','Fontsize', 15, 'Fontweight','b');
ylabel('Temperature coolant outlet [C]','FontName','Arial','Fontsize', 15, 'Fontweight','b');


%Componenti reattivit? 
figure(15)
set(15,'name',['Ramp ' reat ' pcm, Reactivity components']);
hold on;
plot(b_reat_UP(:,1),b_reat_UP(:,4)*1e5,'r','Linewidth',3);
plot(b_reat_UP(:,1),b_reat_UP(:,5)*1e5,'b','Linewidth',3);
plot(b_reat_UP(:,1),b_reat_UP(:,6)*1e5,'m','Linewidth',3);
plot(b_reat_UP(:,1),b_reat_UP(:,7)*1e5,'g','Linewidth',3);
plot(b_reat_UP(:,1),b_reat_UP(:,8)*1e5,'y','Linewidth',3);
hold off;
set(gca,'XLim',[50 1000], 'FontName','Arial','Fontsize', 15, 'Fontweight','b');
grid on;
xlabel('Time [s]','FontName','Arial','Fontsize', 15, 'Fontweight','b');
ylabel('\rho [pcm]','FontName','Arial','Fontsize', 15, 'Fontweight','b');
legend('Doppler','Axial','Coolant','Radial','External','Location','Best');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

