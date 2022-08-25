function [] = plot_V_batch(tExp,yExp,t,y,parameters, title, r2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


if nargin < 7
    r2= zeros(1,4);
    title = 'Veillonella';
elseif  nargin < 6
    title = 'Veillonella';
end
% options 
fs=18;
lw=2;
ms=30;
% get data out of cells
%from simulation 
y1 = y{1};
y2 = y{2};

t1 = t{1};
t2 = t{2};

%from data
yExp1 = yExp{1};
yExp2 = yExp{2};

tExp1 = tExp{1}; 
tExp2 = tExp{2}; 


if length(yExp)>2 
    y3 = y{3};
    t3 = t{3};
    yExp3 = yExp{3};
    tExp3 = tExp{3}; 
end

pos_lac = (strcmp(parameters.compoundAbb,'Slac'));
pos_pro = (strcmp(parameters.compoundAbb,'Spro'));
pos_ac = (strcmp(parameters.compoundAbb,'Sac'));
pos_bm = (strcmp(parameters.compoundAbb,'Xlac'));
pos_ye = (strcmp(parameters.compoundAbb,'Sye'));

figure
set(gcf,'color','w','Position',  [300, 300, 1000, 700]); %background color and size fram of plot
%% data set 1
subplot(2,2,1)
plot(t1,y1(:,pos_lac),'b')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Lactate (g COD/L)')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on
plot(tExp1,yExp1(:,1),'bx')



subplot(2,2,2)
plot(t1,y1(:,pos_pro),'b')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Proprionate (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on
plot(tExp1,yExp1(:,2),'bx')

subplot(2,2,3)
plot(t1,y1(:,pos_ac),'b')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Acetate (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on
plot(tExp1,yExp1(:,3),'bx')

subplot(2,2,4)
plot(t1,y1(:,pos_bm),'b')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('BM Lactate (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on
plot(tExp1,yExp1(:,4),'bx')


%% data Set 2

subplot(2,2,1)
plot(t2,y2(:,pos_lac),'r')
grid on
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Lactate (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')
plot(tExp2,yExp2(:,1),'rx')
text(30,max(y2(:,pos_lac))/4,['r^2=' ,num2str(r2(1))])

subplot(2,2,2)
plot(t2,y2(:,pos_pro),'r')
grid on 
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Proprionate (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
plot(tExp2,yExp2(:,2),'rx')
text(30,max(y2(:,pos_pro))/4,['r^2=' ,num2str(r2(2))])

subplot(2,2,3)
plot(t2,y2(:,pos_ac),'r')
grid on
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Acetate (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
plot(tExp2,yExp2(:,3),'rx')
text(30,max(y2(:,pos_ac))/4,['r^2=' ,num2str(r2(3))])

subplot(2,2,4)
plot(t2,y2(:,pos_bm),'r')
grid on
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('BM Lactate (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
plot(tExp2,yExp2(:,4),'rx')
text(30,max(y2(:,pos_bm))/4, ['r^2=' ,num2str(r2(4))])



%% data Set 3 (optional)
if length(yExp) > 2
subplot(2,2,1)
plot(t3,y3(:,pos_lac),'k')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Glucose (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')
plot(tExp3,yExp3(:,1),'kx')

subplot(2,2,2)
plot(t3,y3(:,pos_pro),'k')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Lactate (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
plot(tExp3,yExp3(:,2),'kx')

subplot(2,2,3)
plot(t3,y3(:,pos_ac),'k')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Yeast extract (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')

subplot(2,2,4)
plot(t3,y3(:,pos_bm),'k')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('BM glucose (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
plot(tExp3,yExp3(:,3),'kx')
end

%sgt = sgtitle(title,'Color','black');
%sgt.FontSize = 20;

%% if you want to plot yeast as well 
% figure 
% plot(t1,y1(:,pos_ye))
% grid on
% xlabel('Time (h)')
% ylabel('Yeast concentration (gCOD/L)')
end

