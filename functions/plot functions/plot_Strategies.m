function [] = plot_Strategies(t,y,parameters)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% options
fs=18;  % font size

% find positions 
pos_glu = (strcmp(parameters.compoundAbb,'Sglu'));
pos_lac = (strcmp(parameters.compoundAbb,'Slac'));
pos_pro = (strcmp(parameters.compoundAbb,'Spro'));



%% get the time of spikes and reed extractions on the differnt graphs
tSpike = parameters.spikeFlowrate.tSpikes-0.08; %spike in sim starts a bit earlier
tReedStart = parameters.Reed.tStart-0.1 ;
tReedEnd = parameters.Reed.tEnd-0.1 ;

matrixTime = repmat(t,1,length(tSpike));
matrixExpTime = repmat(tSpike',length(t),1);
[~,indexExpTimeSpike] = min(abs(matrixExpTime-matrixTime));


 
%% size marker 
mkSize = 8;
mkSizeStar = 16;
figure
set(gcf,'color','w', 'Position',  [300, 300, 900, 700]); %background color and size fram of plot

% set(0,'DefaultLineLineWidth',1,'DefaultAxesFontSize',20,...
%     'DefaultAxesXGrid','off','DefaultAxesYGrid','off',...
%     'DefaultLineMarkerSize',13.5,...
%     'DefaultAxesFontName', 'Calibri', ...  %'Century Schoolbook',...
%     'DefaultAxesXcolor', 0.25*[1,1,1],'DefaultAxesYcolor', 0.25*[1,1,1]);

%% letters of subplots
nIDs = 4;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); % {'(a)','(b)','(c)','(d)'}
letterSize = 18;
positions = [pos_glu, pos_lac, pos_pro]; 
yNames = {'Glucose (gCOD/L)', 'Lactate (gCOD/L)' ,'Proprionate (gCOD/L)', 'Conc in reservoir (gCOD/L)'};

for i = 1:nIDs
    if i < 4
        subplot(2,2,i)
        yElement = y(:,positions(:,i));
        plot(t,yElement,'b')
        xlabel('Time (h)')%,'FontSize',fs,'FontWeight','bold')
        ylabel(yNames{i})
        hold on
        
        %% draw extraction shading
        bd = ylim;
        plot(tSpike(1),yElement(indexExpTimeSpike(1)),'p','MarkerFaceColor','r','MarkerSize' ,mkSizeStar)
        plot(tSpike(2:end),yElement(indexExpTimeSpike(2:end)),'kd','MarkerFaceColor','k','MarkerSize' ,mkSize)
        %a= annotation('arrow',10,5);
        for j = 1:length(tReedStart)
            ptchidx = (t >= tReedStart(j)) & (t <= tReedEnd(j));
            patch([t(ptchidx)' fliplr(t(ptchidx)')], [(zeros(size(yElement(ptchidx)))+bd(2))' zeros(size(yElement(ptchidx)))'], 'g', 'FaceAlpha',0.3, 'EdgeColor','none')
        end
        
        %% text on graph
        x1 = xlim;
        y1 = ylim;
        xt = x1(2)*0.1;
        yt = y1(2)*0.9;
        str = charlbl{i};
        text(xt,yt,str,'FontSize',letterSize,'FontWeight', 'bold','FontName','Calibri')
        grid on
        
    else
        %% extraction
        reservoir = zeros(length(t),14);
        helpReservoir = zeros(1,14);
        volReservoir = 1; %L
        
        for j = 2:length(t)
            states = y(j,:)';
            [~,~,transferReed]=f_mass_balances(t(j),states,parameters);
            helpReservoir = helpReservoir + transferReed'*(t(j)-t(j-1))*states(1)/volReservoir;
            reservoir(j,:) =helpReservoir;
            
        end
        %% plot results
        subplot(2,2,i)
        compoundAbb=parameters.compoundAbb;
        lacR = reservoir(:,strcmp(compoundAbb,'Slac'));
        proR = reservoir(:,strcmp(compoundAbb,'Spro'));
        aceR = reservoir(:,strcmp(compoundAbb,'Sac'));
        
        plot(t,lacR,'b-.',t,proR,'r--',t,aceR,'m:','LineWidth',2)
        legend('Lactate','Propionate','Acetate')
        ylabel(yNames{i})
        xlabel('Time (h)')
        % text on graph
        x1 = xlim;
        y1 = ylim;
        xt = x1(2)*0.1;
        yt = y1(2)*0.9;
        str = charlbl{i};
        text(xt,yt,str,'FontSize',letterSize,'FontWeight', 'bold','FontName','Calibri')
        grid on
        
    end
end

hold off

end

