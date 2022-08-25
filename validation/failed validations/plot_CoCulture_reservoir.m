function [] = plot_CoCulture_reservoir(tExp,yExp,t,y,parameters,r2,ReservoirData)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



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
set(gcf,'color','w', 'Position',  [300, 300, 1000, 700]); %background color

set(0,'DefaultLineLineWidth',1,'DefaultAxesFontSize',18,...
    'DefaultAxesXGrid','off','DefaultAxesYGrid','off',...
    'DefaultLineMarkerSize',13.5,...
    'DefaultAxesFontName', 'Calibri') ...  %'Century Schoolbook',...
    %'DefaultAxesXcolor', 0.25*[1,1,1],'DefaultAxesYcolor', 0.25*[1,1,1]);

%% letters of subplots
nIDs = 6;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); % {'(a)','(b)','(c)','(d)'}
text(0.025,0.95,charlbl{1},'Units','normalized','FontSize',12)
letterSize = 20;

%% Glucose subplot
subplot(2,3,1)
yGlu = y(:,pos_glu);
plot(t,y(:,pos_glu),'b')
xlabel('Time (h)')%,'FontSize',fs,'FontWeight','bold')
ylabel('Glucose (gCOD/L)')
%set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on
plot(tExp,yExp(:,1),'bx')
% plot(tSpike,y(indexExpTimeSpike,pos_glu),'kd','MarkerFaceColor','k','MarkerSize' ,mkSize)
% plot(tReedStart,y(indexExpTimeReedStart,pos_glu),'kd','MarkerFaceColor','g','MarkerSize' ,mkSize)
% plot(tReedEnd,y(indexExpTimeReedEnd,pos_glu),'kd','MarkerFaceColor','r','MarkerSize' ,mkSize)

%% draw extraction 
bd = ylim;
plot(tSpike,yGlu(indexExpTimeSpike),'kd','MarkerFaceColor','k','MarkerSize' ,mkSize)
for i= 1:length(tReedStart)
    ptchidx = (t >= tReedStart(i)) & (t <= tReedEnd(i));
    patch([t(ptchidx)' fliplr(t(ptchidx)')], [(zeros(size(yGlu(ptchidx)))+bd(2))' zeros(size(yGlu(ptchidx)))'], 'g', 'FaceAlpha',0.3, 'EdgeColor','none') 
end
% text on graph
x1 = xlim;
y1 = ylim;
xt = x1(2)*0.1;
yt = y1(2)*0.9;
%str = {charlbl{1},['R^2=' ,num2str(r2(1))]};
str = charlbl{1};
text(xt,yt,str,'FontSize',letterSize,'FontWeight', 'bold','FontName','Calibri')

grid on
%% Lactate subplot
subplot(2,3,2)
yLac = y(:,pos_lac);
plot(t,yLac,'b')
xlabel('Time (h)')%,'FontSize',fs,'FontWeight','bold')
ylabel('Lactate (gCOD/L)')
%set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on
plot(tExp,yExp(:,2),'bx')

%% Extraction 
bd = ylim;
plot(tSpike,yLac(indexExpTimeSpike),'kd','MarkerFaceColor','k','MarkerSize' ,mkSize)
for i= 1:length(tReedStart)
    ptchidx = (t >= tReedStart(i)) & (t <= tReedEnd(i));
    patch([t(ptchidx)' fliplr(t(ptchidx)')], [(zeros(size(yLac(ptchidx)))+bd(2))' zeros(size(yLac(ptchidx)))'], 'g', 'FaceAlpha',0.3, 'EdgeColor','none') 
end

% text on graph
x1 = xlim;
y1 = ylim;
xt = x1(2)*0.1;
yt = y1(2)*0.9;
%str = {charlbl{2},['R^2=' ,num2str(r2(2))]};
str = charlbl{2};
text(xt,yt,str,'FontSize',letterSize,'FontWeight', 'bold','FontName','Calibri')
grid on
%% Proprionate subplot
subplot(2,3,3)
yProp  = y(:,pos_pro);
plot(t,yProp,'b')
xlabel('Time (h)')%,'FontSize',fs,'FontWeight','bold')
ylabel('Proprionate (gCOD/L)')%,'FontSize',fs,'FontWeight','bold')
%set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on
plot(tExp,yExp(:,3),'bx')

%% draw extraction 
bd = ylim;
plot(tSpike,yProp(indexExpTimeSpike),'kd','MarkerFaceColor','k','MarkerSize' ,mkSize)
for i= 1:length(tReedStart)
    ptchidx = (t >= tReedStart(i)) & (t <= tReedEnd(i));
    patch([t(ptchidx)' fliplr(t(ptchidx)')], [(zeros(size(yProp(ptchidx)))+bd(2))' zeros(size(yProp(ptchidx)))'], 'g', 'FaceAlpha',0.3, 'EdgeColor','none') 
end

%% text on graph
x1 = xlim;
y1 = ylim;
xt = x1(2)*0.1;
yt = y1(2)*0.9;
%str = {charlbl{3},['R^2=' ,num2str(r2(3))]};
str = charlbl{3};
text(xt,yt,str,'FontSize',letterSize,'FontWeight', 'bold','FontName','Calibri')
grid on

%% Acetate subplot
subplot(2,3,4)
yAce= y(:,pos_ac);
plot(t,y(:,pos_ac),'b')
xlabel('Time (h)')%,'FontSize',fs,'FontWeight','bold')
ylabel('Acetate (gCOD/L)')%,'FontSize',fs,'FontWeight','bold')
%set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on
plot(tExp,yExp(:,4),'bx')
grid on

%% draw extraction 
bd = ylim;
plot(tSpike,yAce(indexExpTimeSpike),'kd','MarkerFaceColor','k','MarkerSize' ,mkSize)
for i= 1:length(tReedStart)
    ptchidx = (t >= tReedStart(i)) & (t <= tReedEnd(i));
    patch([t(ptchidx)' fliplr(t(ptchidx)')], [(zeros(size(yAce(ptchidx)))+bd(2))' zeros(size(yAce(ptchidx)))'], 'g', 'FaceAlpha',0.3, 'EdgeColor','none') 
end

% text on graph
x1 = xlim;
y1 = ylim;
xt = x1(2)*0.1;
yt = y1(2)*0.9;
%str = {charlbl{4},['R^2=' ,num2str(r2(4))]};
str = charlbl{4};
text(xt,yt,str,'FontSize',letterSize,'FontWeight', 'bold','FontName','Calibri')
%% BM total subplot
subplot(2,3,5)
bm =  y(:,pos_bm_Su)+y(:,pos_bm_lac) ;
plot(t,bm,'b')
xlabel('Time (h)')%,'FontSize',fs,'FontWeight','bold')
ylabel('Total biomass (gCOD/L)')%,'FontSize',fs,'FontWeight','bold')
%set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on
plot(tExp,yExp(:,5),'bx')

% text on graph
x1 = xlim;
y1 = ylim;
xt = x1(2)*0.1;
yt = y1(2)*0.9;
%str = {charlbl{5},['R^2=' ,num2str(r2(5))]};
str = charlbl{5};
text(xt,yt,str,'FontSize',letterSize,'FontWeight', 'bold','FontName','Calibri')
% grid
grid on
%% plot reservoir
subplot(2,3,6) 
load(ReservoirData)
compoundAbb = parameters.compoundAbb;

lacR = reservoir(:,strcmp(compoundAbb,'Slac'));
proR = reservoir(:,strcmp(compoundAbb,'Spro'));
aceR = reservoir(:,strcmp(compoundAbb,'Sac'));

%plot(t,lacR,'b',t,proR,'r',t,aceR,'m')
plot(t,lacR,'b-.',t,proR,'r--',t,aceR,'m:','LineWidth',1.5)
hold on 
plot(tRV,lacRV,'bx',tRV,proRV,'rx',tRV,aceRV,'mx')
  
%lgd = legend('Lactate','Propionate','Acetate');
%lgd.FontSize = 13;
ylabel('Conc in reservoir (gCOD/L)') 
xlabel('Time (h)')
% text on graph
x1 = xlim;
y1 = ylim;
xt = x1(2)*0.1;
yt = y1(2)*0.9;
str = charlbl{6};
text(xt,yt,str,'FontSize',letterSize,'FontWeight', 'bold','FontName','Calibri')
grid on


%%
hold off 
end

