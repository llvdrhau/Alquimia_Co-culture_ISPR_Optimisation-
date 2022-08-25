function [] = f_plot(tExp,yExp,t,y,parameters)

% lt = length(t);
% yExp1 = yExp(:,1:3);
% y1 = y(1:lt,:);
% yExp2 = yExp(:,4:6);
% y2 = y(lt+1:2*lt,:);

% get data out of cells
%from simulation 
y1 = y{1};
y2 = y{2};
y3 = y{3};

t1 = t{1};
t2 = t{2};
t3 = t{3};

%from data
yExp1 = yExp{1};
yExp2 = yExp{2};
yExp3 = yExp{3};

tExp1 = tExp{1}; 
tExp2 = tExp{2}; 
tExp3 = tExp{3}; 

% wat is this? 
rgb=[231 172 23; 77 127 97; 150 150 150]/255;
rgb2=[186 138 18; 46 76 58; 85 85 85]/255;
fs=18;
lw=2;
ms=30;

pos_glu = (strcmp(parameters.compoundAbb,'Sglu'));
pos_lac = (strcmp(parameters.compoundAbb,'Slac'));
pos_bm = (strcmp(parameters.compoundAbb,'Xsu'));
pos_ye = (strcmp(parameters.compoundAbb,'Sye'));
pos_tan = (strcmp(parameters.compoundAbb,'SIn'));

figure
%% Set 1
subplot(3,2,1)
plot(t1,y1(:,pos_glu),'color',rgb2(1,:),'LineWidth',lw)
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Glucose (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'FontWeight','bold')
hold on
plot(tExp1,yExp1(:,1),'.','color',rgb(1,:),'MarkerSize',ms)

subplot(3,2,2)
plot(t1,y1(:,pos_lac),'color',rgb2(1,:),'LineWidth',lw)
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Lactate (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on
plot(tExp1,yExp1(:,2),'.','color',rgb(1,:),'MarkerSize',ms)

subplot(3,2,3)
plot(t1,y1(:,pos_ye),'color',rgb2(1,:),'LineWidth',lw)
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Yeast extract (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on

subplot(3,2,4)
plot(t1,y1(:,pos_bm),'color',rgb2(1,:),'LineWidth',lw)
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('BM glucose (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on
plot(tExp1,yExp1(:,3),'.','color',rgb(1,:),'MarkerSize',ms)

subplot(3,2,5)
plot(t1,y1(:,pos_tan),'color',rgb2(1,:),'LineWidth',lw)
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('TAN (mol N/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
hold on

%% Set 2

subplot(3,2,1)
plot(t2,y2(:,pos_glu),'color',rgb2(2,:),'LineWidth',lw)
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Glucose (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'FontWeight','bold')
plot(tExp2,yExp2(:,1),'.','color',rgb(2,:),'MarkerSize',ms)

subplot(3,2,2)
plot(t2,y2(:,pos_lac),'color',rgb2(2,:),'LineWidth',lw)
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Lactate (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
plot(tExp2,yExp2(:,2),'.','color',rgb(2,:),'MarkerSize',ms)

subplot(3,2,3)
plot(t2,y2(:,pos_ye),'color',rgb2(2,:),'LineWidth',lw)
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Yeast extract (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')

subplot(3,2,4)
plot(t2,y2(:,pos_bm),'color',rgb2(2,:),'LineWidth',lw)
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('BM glucose (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
plot(tExp2,yExp2(:,3),'.','color',rgb(2,:),'MarkerSize',ms)

subplot(3,2,5)
plot(t2,y2(:,pos_tan),'color',rgb2(2,:),'LineWidth',lw)
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('TAN (mol N/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')


%% Set 3

subplot(3,2,1)
plot(t3,y3(:,pos_glu),'color',rgb2(3,:),'LineWidth',lw)
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Glucose (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'FontWeight','bold')
plot(tExp3,yExp3(:,1),'.','color',rgb(3,:),'MarkerSize',ms)

subplot(3,2,2)
plot(t3,y3(:,pos_lac),'color',rgb2(3,:),'LineWidth',lw)
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Lactate (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
plot(tExp3,yExp3(:,2),'.','color',rgb(3,:),'MarkerSize',ms)

subplot(3,2,3)
plot(t3,y3(:,pos_ye),'color',rgb2(3,:),'LineWidth',lw)
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('Yeast extract (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')

subplot(3,2,4)
plot(t3,y3(:,pos_bm),'color',rgb2(3,:),'LineWidth',lw)
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('BM glucose (g COD/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')
plot(tExp3,yExp3(:,3),'.','color',rgb(3,:),'MarkerSize',ms)

subplot(3,2,5)
plot(t3,y3(:,pos_tan),'color',rgb2(3,:),'LineWidth',lw)
xlabel('Time (h)','FontSize',fs,'FontWeight','bold')
ylabel('TAN (mol N/L)','FontSize',fs,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',12,'FontWeight','bold')


end

