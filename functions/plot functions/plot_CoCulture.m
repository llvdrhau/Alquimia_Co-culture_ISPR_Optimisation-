function [] = plot_CoCulture(tExp,yExp,t,y,parameters, title, r2, plotYeast)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


if nargin < 6
    title = 'Co-Culture';
    r2 = zeros(1,5);
    plotYeast = 0;
    
elseif nargin < 7
    r2 = zeros(1,5);
    plotYeast = 0;
    
elseif nargin < 8
    plotYeast = 0;
  
end

% options
fs=18;  % font size


pos_glu = (strcmp(parameters.compoundAbb,'Sglu'));
pos_lac = (strcmp(parameters.compoundAbb,'Slac'));
pos_pro = (strcmp(parameters.compoundAbb,'Spro'));
pos_ac = (strcmp(parameters.compoundAbb,'Sac'));
pos_bm_lac = (strcmp(parameters.compoundAbb,'Xlac'));
pos_bm_Su = (strcmp(parameters.compoundAbb,'Xsu'));
pos_ye = (strcmp(parameters.compoundAbb,'Sye'));



%% get the time of spikes and reed extractions on the differnt graphs
tSpike = parameters.spikeFlowrate.tSpikes-0.2; %spike in sim starts a bit earlier 
tReedStart = parameters.Reed.tStart ; 
tReedEnd = parameters.Reed.tEnd ;  


matrixTime = repmat(t,1,length(tSpike));
matrixExpTime = repmat(tSpike',length(t),1);
[~,indexExpTimeSpike] = min(abs(matrixExpTime-matrixTime));


matrixTime = repmat(t,1,length(tReedStart));
matrixExpTime = repmat(tReedStart',length(t),1);
[~,indexExpTimeReedStart] = min(abs(matrixExpTime-matrixTime));


matrixTime = repmat(t,1,length(tReedEnd));
matrixExpTime = repmat(tReedEnd',length(t),1);
[~,indexExpTimeReedEnd] = min(abs(matrixExpTime-matrixTime));

%% marker size
mkSize = 8;
set(gcf,'color','w'); %background color
%% Glucose subplot
subplot(2,3,1)
plot(t,y(:,pos_glu),'b')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Glucose (gCOD/L)')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on
plot(tExp,yExp(:,1),'bx')
plot(tSpike,y(indexExpTimeSpike,pos_glu),'kd','MarkerFaceColor','k','MarkerSize' ,mkSize)
plot(tReedStart,y(indexExpTimeReedStart,pos_glu),'kd','MarkerFaceColor','g','MarkerSize' ,mkSize)
plot(tReedEnd,y(indexExpTimeReedEnd,pos_glu),'kd','MarkerFaceColor','r','MarkerSize' ,mkSize)

text(max(t)*0.1,max(y(:,pos_glu))*1.2, ['r^2=' ,num2str(r2(1))],'FontSize',12)
grid on
%% Lactate subplot
subplot(2,3,2)
plot(t,y(:,pos_lac),'b')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Lactate (gCOD/L)')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on
plot(tReedStart,y(indexExpTimeReedStart,pos_lac),'kd','MarkerFaceColor','g','MarkerSize' ,mkSize)
plot(tReedEnd,y(indexExpTimeReedEnd,pos_lac),'kd','MarkerFaceColor','r','MarkerSize' ,mkSize)
plot(tExp,yExp(:,2),'bx')
text(max(t)*0.1,max(y(:,pos_lac))*0.9, ['r^2=' ,num2str(r2(2))],'FontSize',12)
grid on
%% Proprionate subplot
subplot(2,3,3)
plot(t,y(:,pos_pro),'b')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Proprionate (gCOD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on
plot(tExp,yExp(:,3),'bx')
plot(tReedStart,y(indexExpTimeReedStart,pos_pro),'kd','MarkerFaceColor','g','MarkerSize' ,mkSize)
plot(tReedEnd,y(indexExpTimeReedEnd,pos_pro),'kd','MarkerFaceColor','r','MarkerSize' ,mkSize)
text(max(t)*0.1,max(y(:,pos_pro))*0.9, ['r^2=' ,num2str(r2(3))],'FontSize',12)
grid on

%% Acetate subplot
subplot(2,3,4)
plot(t,y(:,pos_ac),'b')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Acetate (gCOD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on
plot(tExp,yExp(:,4),'bx')
plot(tReedStart,y(indexExpTimeReedStart,pos_ac),'kd','MarkerFaceColor','g','MarkerSize' ,mkSize)
plot(tReedEnd,y(indexExpTimeReedEnd,pos_ac),'kd','MarkerFaceColor','r','MarkerSize' ,mkSize)
text(max(t)*0.1,max(y(:,pos_ac))*0.9, ['r^2=' ,num2str(r2(4))],'FontSize',12)
grid on

%% BM total subplot
subplot(2,3,5)
bm =  y(:,pos_bm_Su)+y(:,pos_bm_lac) ;
plot(t,bm,'b')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Total biomass (gCOD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on
plot(tExp,yExp(:,5),'bx')
text(max(t)*0.1,max(bm)*0.9, ['r^2=' ,num2str(r2(5))],'FontSize',12)
grid on

%% BMLac subplot
subplot(2,3,6)
bm = y(:,pos_bm_lac) ;
plot(t,bm,'b')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Biomass Veillionella (gCOD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
% hold on
% plot(tExp,yExp(:,5),'bx')
grid on

%%
%% transfer to Reservoir
% reservoir = zeros(length(t),14); 
% volReservoir = 1; %L 
% for i = 2:length(t)
%     states = y(i,:)';
%     [~,~,transferReed]=f_mass_balances(t(i),states,parameters); 
%     reservoir(i,:) = reservoir(i-1,:) + transferReed'*(t(i)-t(i-1))*states(1)/volReservoir;
% end
%  
% lacR = reservoir(:,strcmp(compoundAbb,'Slac'));
% proR = reservoir(:,strcmp(compoundAbb,'Spro'));
% aceR = reservoir(:,strcmp(compoundAbb,'Sac'));
% figure
% plot(t,lacR,t,proR,t,aceR)
% hold on 
% tRV=data1(:,strcmp(Names,'timeRV'));
% gluRV = data1(:,strcmp(Names,'GluRV'));
% lacRV = data1(:,strcmp(Names,'LacRV'));
% proRV = data1(:,strcmp(Names,'ProRV'));
% acRV = data1(:,strcmp(Names,'AceRV'));
% plot(tRV,lacRV,'bx',tRV,proRV,'rx',tRV,acRV,'yx')
% 
% legend('Lac','Pro','Ace')
% ylabel('concentration g/L') 
% xlabel('Time h')

%% title
sgt = sgtitle(title,'Color','black');
sgt.FontSize = 20;

%% YE plot
if plotYeast 
plot(t,y(:,pos_ye))
grid on
xlabel('Time (h)')
ylabel('Yeast concentration (gCOD/L)')
end


%% vol
figure 
plot(t,y(:,1))
xlabel('Time (h)')
ylabel('Volume (L)')
%%
hold off 
end

