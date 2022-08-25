function [] = plot_CoCulture_simulation(t,y,parameters, title ,r2, plotYeast)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


if nargin < 4
    title = 'Co-Culture';
    r2 = zeros(1,5);
    plotYeast = 0;
elseif nargin < 5
    r2 = zeros(1,5);
    plotYeast = 0;
elseif nargin < 6
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
tSpike = parameters.spikeFlowrate.tSpikes-0.08; %spike in sim starts a bit earlier
tReedStart = parameters.Reed.tStart-0.1 ;
tReedEnd = parameters.Reed.tEnd-0.1 ;

%t = union(t,tSpike);
%indexCell = cell(nDataExp,nExp);   % preallocation, use if more than  one
%dataset

matrixTime = repmat(t,1,length(tSpike));
matrixExpTime = repmat(tSpike',length(t),1);
[~,indexExpTimeSpike] = min(abs(matrixExpTime-matrixTime));

%t = union(t,tReedStart);
%indexCell = cell(nDataExp,nExp);   % preallocation, use if more than  one
%dataset

% matrixTime = repmat(t,1,length(tReedStart));
% matrixExpTime = repmat(tReedStart',length(t),1);
% [~,indexExpTimeReedStart] = min(abs(matrixExpTime-matrixTime));
% 
% %t = union(t,tReedEnd);
% %indexCell = cell(nDataExp,nExp);   % preallocation, use if more than  one
% %dataset
% 
% matrixTime = repmat(t,1,length(tReedEnd));
% matrixExpTime = repmat(tReedEnd',length(t),1);
% [~,indexExpTimeReedEnd] = min(abs(matrixExpTime-matrixTime));
 
%% size marker 
mkSize = 8;
lineWidth = 4;
set(gcf,'color','w', 'Position',  [300, 300, 1000, 700]); %background color
%% letters of subplots
nIDs = 6;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); % {'(a)','(b)','(c)','(d)'}

%% Glucose subplot
subplot(2,3,1)
yGlu = y(:,pos_glu);
plot(t,yGlu,'b')
xlabel('Time (h)')%,'FontSize',fs,'FontWeight','bold')
ylabel('Glucose (gCOD/L)')
hold on 

%% draw extraction 
bd = ylim;
plot(tSpike(1),yGlu(indexExpTimeSpike(1)),'p','MarkerFaceColor','r','MarkerSize' ,mkSize)
plot(tSpike(2:end),yGlu(indexExpTimeSpike(2:end)),'kd','MarkerFaceColor','k','MarkerSize' ,mkSize)
for i= 1:length(tReedStart)
    ptchidx = (t >= tReedStart(i)) & (t <= tReedEnd(i));
    patch([t(ptchidx)' fliplr(t(ptchidx)')], [(zeros(size(yGlu(ptchidx)))+bd(2))' zeros(size(yGlu(ptchidx)))'], 'g', 'FaceAlpha',0.3, 'EdgeColor','none') 
end

%% text on graph
x1 = xlim;
y1 = ylim;
xt = x1(2)*0.1;
yt = y1(2)*0.9;
str = charlbl{1};
text(xt,yt,str,'FontSize',12)
grid on

%% Lactate subplot
subplot(2,3,2)
yLac = y(:,pos_lac);
plot(t,yLac,'b')
xlabel('Time (h)')%,'FontSize',fs,'FontWeight','bold')
ylabel('Lactate (gCOD/L)')
hold on 
grid on
%% Extraction 
bd = ylim;
plot(tSpike(1),yLac(indexExpTimeSpike(1)),'p','MarkerFaceColor','r','MarkerSize' ,mkSize)
plot(tSpike(2:end),yLac(indexExpTimeSpike(2:end)),'kd','MarkerFaceColor','k','MarkerSize' ,mkSize)
for i= 1:length(tReedStart)
    ptchidx = (t >= tReedStart(i)) & (t <= tReedEnd(i));
    patch([t(ptchidx)' fliplr(t(ptchidx)')], [(zeros(size(yLac(ptchidx)))+bd(2))' zeros(size(yLac(ptchidx)))'], 'g', 'FaceAlpha',0.3, 'EdgeColor','none') 
end

%% text on graph
x1 = xlim;
y1 = ylim;
xt = x1(2)*0.1;
yt = y1(2)*0.9;
str = charlbl{2};
text(xt,yt,str,'FontSize',12)
grid on

%% Proprionate subplot
subplot(2,3,3)
yProp = y(:,pos_pro);
plot(t,yProp,'b')
xlabel('Time (h)')%,'FontSize',fs,'FontWeight','bold')
ylabel('Proprionate (gCOD/L)')%,'FontSize',fs,'FontWeight','bold')
hold on, 

%% draw extraction 
bd = ylim;
plot(tSpike(1),yProp(indexExpTimeSpike(1)),'p','MarkerFaceColor','r','MarkerSize' ,mkSize)
plot(tSpike(2:end),yProp(indexExpTimeSpike(2:end)),'kd','MarkerFaceColor','k','MarkerSize' ,mkSize)
for i= 1:length(tReedStart)
    ptchidx = (t >= tReedStart(i)) & (t <= tReedEnd(i));
    patch([t(ptchidx)' fliplr(t(ptchidx)')], [(zeros(size(yProp(ptchidx)))+bd(2))' zeros(size(yProp(ptchidx)))'], 'g', 'FaceAlpha',0.3, 'EdgeColor','none') 
end

%% text on graph
x1 = xlim;
y1 = ylim;
xt = x1(2)*0.1;
yt = y1(2)*0.9;
str = {charlbl{3}};
text(xt,yt,str,'FontSize',12)
grid on

grid on

%% Acetate subplot
subplot(2,3,4)
yAce = y(:,pos_ac);
plot(t,yAce,'b')
xlabel('Time (h)')%,'FontSize',fs,'FontWeight','bold')
ylabel('Acetate (gCOD/L)')%,'FontSize',fs,'FontWeight','bold')
hold on, 
grid on
%% text on graph
x1 = xlim;
y1 = ylim;
xt = x1(2)*0.1;
yt = y1(2)*0.9;
str = {charlbl{4}};
text(xt,yt,str,'FontSize',12)
grid on

%% draw extraction 
bd = ylim;
plot(tSpike(1),yAce(indexExpTimeSpike(1)),'p','MarkerFaceColor','r','MarkerSize' ,mkSize)
plot(tSpike(2:end),yAce(indexExpTimeSpike(2:end)),'kd','MarkerFaceColor','k','MarkerSize' ,mkSize)
for i= 1:length(tReedStart)
    ptchidx = (t >= tReedStart(i)) & (t <= tReedEnd(i));
    patch([t(ptchidx)' fliplr(t(ptchidx)')], [(zeros(size(yAce(ptchidx)))+ bd(2))' zeros(size(yAce(ptchidx)))'], 'g', 'FaceAlpha',0.3, 'EdgeColor','none') 
end


%% BMSu subplot
subplot(2,3,5)
bm =  y(:,pos_bm_Su)+y(:,pos_bm_lac);
plot(t,bm,'b')
xlabel('Time (h)')%,'FontSize',fs,'FontWeight','bold')
ylabel('Total Biomass (gCOD/L)')%,'FontSize',fs,'FontWeight','bold')
%set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
grid on
% text on graph
x1 = xlim;
y1 = ylim;
xt = x1(2)*0.1;
yt = y1(2)*0.9;
str = {charlbl{5}};
text(xt,yt,str,'FontSize',12)
grid on


%% BMLac subplot
% subplot(2,3,6)
% bm = y(:,pos_bm_lac) ;
% plot(t,bm,'b')
% xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
% ylabel('Biomass Veillionella (gCOD/L)','FontSize',fs,'FontWeight','bold')
% set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
% grid on

%% title
% sgt = sgtitle(title,'Color','black');
% sgt.FontSize = 20;


%% reservoir 
%% transfer to Reservoir old code 
% reservoir = zeros(length(t),14); 
% volReservoir = 1; %L 
% for i = 2:length(t)
%     time = t(i);
%     states = y(i,:)';
%     [~,~,transferReed]=f_mass_balances(t(i),states,parameters); 
%     reservoir(i,:) = reservoir(i-1,:) + transferReed'*(t(i)-t(i-1))*states(1)/volReservoir;
% end

% compoundAbb=parameters.compoundAbb; 
% lacR = reservoir(:,strcmp(compoundAbb,'Slac'));
% proR = reservoir(:,strcmp(compoundAbb,'Spro'));
% aceR = reservoir(:,strcmp(compoundAbb,'Sac'));
% figure
% plot(t,lacR,t,proR,t,aceR)
% 
% legend('Lac','Pro','Ace')
% ylabel('concentration g/L') 
% xlabel('Time h')
% hold off 

%% alternative methode of caluculation reservoir 
% 
% % reshape data so time steps taken are the same... namely 0.1 
% timeNew = (0:0.1:t(end))';
% [~,nStates] = size(y);
% newState = zeros(length(timeNew),nStates);  % preallocation 
% 
% for i = 1:nStates
% newState(:,i) = interp1(t,y(:,i),timeNew);
% end 
% 
% % calcualte reservoir
% reservoir = zeros(length(timeNew),nStates); 
% helpReservoir = zeros(1,nStates); 
% volReservoir = 1; %L 
% 
% for i = 2:length(timeNew)    
%     states = newState(i,:)';
%     [~,~,transferReed]=f_mass_balances(timeNew(i),states,parameters); 
%     helpReservoir = helpReservoir + transferReed'*(timeNew(i)-timeNew(i-1))*states(1)/volReservoir;
%     reservoir(i,:) =helpReservoir;    
% end

%% extraction test 2
reservoir = zeros(length(t),14); 
helpReservoir = zeros(1,14); 
volReservoir = 1; %L 

for i = 2:length(t)    
    states = y(i,:)';
    [~,~,transferReed]=f_mass_balances(t(i),states,parameters); 
    helpReservoir = helpReservoir + transferReed'*(t(i)-t(i-1))*states(1)/volReservoir;
    reservoir(i,:) =helpReservoir;
    
end
%% plot results 
subplot(2,3,6)
compoundAbb=parameters.compoundAbb; 
lacR = reservoir(:,strcmp(compoundAbb,'Slac'));
proR = reservoir(:,strcmp(compoundAbb,'Spro'));
aceR = reservoir(:,strcmp(compoundAbb,'Sac'));

%set(gcf,'color','w'); %background color
%plot(timeNew,lacR,timeNew,proR,timeNew,aceR)
plot(t,lacR,t,proR,t,aceR)

legend('Lactate','Propionate','Acetate')
ylabel('Concentration in reservoir (gCOD/L)') 
xlabel('Time (h)')
% text on graph
x1 = xlim;
y1 = ylim;
xt = x1(2)*0.1;
yt = y1(2)*0.9;
str = {charlbl{6}};
text(xt,yt,str,'FontSize',12)
grid on

hold off 


%% YE plot
if plotYeast 
plot(t,y(:,pos_ye))
grid on
xlabel('Time (h)')
ylabel('Yeast concentration (gCOD/L)')
end


%% vol
% figure 
% plot(t,y(:,1))
% xlabel('Time (h)')
% ylabel('Volume (L)')
end

