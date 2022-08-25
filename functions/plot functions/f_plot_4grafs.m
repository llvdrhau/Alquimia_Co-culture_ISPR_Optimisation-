function [] = f_plot_4grafs(tExp,yExp,t,y,parameters)

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

pos_glu = (strcmp(parameters.compoundAbb,'Sglu'));
pos_lac = (strcmp(parameters.compoundAbb,'Slac'));
pos_bm = (strcmp(parameters.compoundAbb,'Xsu'));
pos_ye = (strcmp(parameters.compoundAbb,'Sye'));

figure
%% Set 1
subplot(2,2,1)
plot(t1,y1(:,pos_glu),'b')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Glucose (g COD/L)')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on
plot(tExp1,yExp1(:,1),'bx')

subplot(2,2,2)
plot(t1,y1(:,pos_lac),'b')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Lactate (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on
plot(tExp1,yExp1(:,2),'bx')

subplot(2,2,3)
plot(t1,y1(:,pos_ye),'b')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Yeast extract (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on

subplot(2,2,4)
plot(t1,y1(:,pos_bm),'b')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('BM glucose (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on
plot(tExp1,yExp1(:,3),'bx')


%% Set 2

subplot(2,2,1)
plot(t2,y2(:,pos_glu),'g')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Glucose (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'FontWeight','bold')
plot(tExp2,yExp2(:,1),'gx')

subplot(2,2,2)
plot(t2,y2(:,pos_lac),'g')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Lactate (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
plot(tExp2,yExp2(:,2),'gx')
subplot(2,2,3)
plot(t2,y2(:,pos_ye),'g')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Yeast extract (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')

subplot(2,2,4)
plot(t2,y2(:,pos_bm),'g')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('BM glucose (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
plot(tExp2,yExp2(:,3),'gx')

%% Set 3
if length(yExp) > 2
subplot(2,2,1)
plot(t3,y3(:,pos_glu),'k')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Glucose (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'FontWeight','bold')
plot(tExp3,yExp3(:,1),'kx')

subplot(2,2,2)
plot(t3,y3(:,pos_lac),'k')
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Lactate (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
plot(tExp3,yExp3(:,2),'kx')

subplot(2,2,3)
plot(t3,y3(:,pos_ye),'k')
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

end

