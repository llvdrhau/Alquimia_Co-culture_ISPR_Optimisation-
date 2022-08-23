function [] = f_plot_montecarlo_Veillionella(Matrix1,Matrix2,tExp,yExp)
%% This script plots the solution space determined in the Monte Carlo method
% Alberte Regueira
% Department of Chemical Engineering, USC
% November 2019

t=linspace(0,tExp(end),50);
t = union(t,tExp);


matrix1_glu = Matrix1.glu;
matrix1_lac = Matrix1.lac;
matrix1_pro = Matrix1.pro;
matrix1_ace = Matrix1.ace;
matrix1_ye = Matrix1.ye;
matrix1_bm = Matrix1.bm_lac;
matrix1_tan = Matrix1.tan;

matrix2_glu = Matrix2.glu;
matrix2_lac = Matrix2.lac;
matrix2_pro = Matrix2.pro;
matrix2_ace = Matrix2.ace;
matrix2_ye = Matrix2.ye;
matrix2_bm = Matrix2.bm_lac;
matrix2_tan = Matrix2.tan;


compounds = {'lac','pro','ace','bm','ye'}; % add yeat to the plots? 
%compounds = {'lac','pro','ace','bm','ye'};
y_axis = {'Lactate (gCOD/L)','Proprionate (gCOD/L)','Acetate (gCOD/L)','Biomass (gCOD/L)','Yeast extract (gCOD/L)'};
%graph_titles = {'Lactate','proprionate','Acetate','Biomass','Yeast extract','TAN'};

% propertie of graph 
set(gcf,'color','w', 'Position',  [300, 300, 1000, 700]); %background color
set(0,'DefaultLineLineWidth',1,'DefaultAxesFontSize',20,...
    'DefaultAxesXGrid','off','DefaultAxesYGrid','off',...
    'DefaultLineMarkerSize',13.5,...
    'DefaultAxesFontName','Calibri',...
    'DefaultAxesXcolor', 0.25*[1,1,1],'DefaultAxesYcolor', 0.25*[1,1,1]);

rgb = [186 138 18; 46 76 58; 132 132 132]/255;
rgb_area = [255 220 122; 157 212 179; 179 179 179]/255;

%% letters of subplots
nIDs = 5;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); % {'(a)','(b)','(c)','(d)'}

%% plots 
j= 2;
for i=1:length(compounds)-1
    
    %subplot(3,2,i)
    subplot(2,2,i)
    
    eval(horzcat('matrix=matrix',num2str(j),'_',compounds{i},';'))
    
    UP = prctile(matrix,97.5,2);
    DOWN =prctile(matrix,2.5,2);
    y = mean(matrix,2);
    
    p = plot(t,y,t,UP,t,DOWN,'g');
    YLIM = get(gca,'YLim');
    delete(p);
    %YLIM(1) = 0;
    a1 = area(t,UP,min(YLIM));
    hold on;
    set(a1,'LineStyle','none');     set(a1,'FaceColor',rgb_area(j,:));
    a2 = area(t,DOWN,min(YLIM));
    set(a2,'LineStyle','none');     set(a2,'FaceColor',[1 1 1]);
    
    hold on
    plot(t,y,'color',rgb(j,:))
    
    if i<5
        plot(tExp,yExp(:,i),'d','color',[0 0.4470 0.7410],'MarkerFaceColor',[0 0.4470 0.7410],'MarkerSize',8) %rgb(1,:)
        hold on 
        plot(tExp,yExp(:,4+i),'.','color',rgb(2,:),'MarkerSize',25)
        %legend('Monte Carlo interval','','Simulation','dataset 2.1','dataset 2.2')
    end
    
    ylabel(y_axis(i))
    xlabel('Time (h)')
    xlim([t(1) t(end)])
    %title(graph_titles{i})
    set(gcf,'color','w'); %background color
    set(gca,'Layer','top','XGrid','on','YGrid','on');
   
    % text on graph
    x1 = xlim;
    y1 = ylim;
    xt = x1(2)*0.1;
    yt = y1(2)*0.8;
    str = charlbl{i};
    text(xt,yt,str,'FontSize',20, 'FontWeight', 'bold','FontName','Calibri')
end



% for j=1:2
%     figure
%     if j==1
%         yEx = yExp(:,1:4); %if you want to seperate them 
%     elseif j==2
%         yEx = yExp(:,5:end);
%     end
%     for i=1:length(compounds)
%         
%         subplot(3,2,i)
%         
%         eval(horzcat('matrix=matrix',num2str(j),'_',compounds{i},';'))
%         
%         UP = prctile(matrix,97.5,2);
%         DOWN =prctile(matrix,2.5,2);
%         y = mean(matrix,2);
%         
%         p = plot(t,y,t,UP,t,DOWN,'g');
%         YLIM = get(gca,'YLim');
%         delete(p);
%         %YLIM(1) = 0;
%         a1 = area(t,UP,min(YLIM));
%         hold on;
%         set(a1,'LineStyle','none');     set(a1,'FaceColor',rgb_area(j,:));
%         a2 = area(t,DOWN,min(YLIM));
%         set(a2,'LineStyle','none');     set(a2,'FaceColor',[1 1 1]);
%         
%         hold on
%         plot(t,y,'color',rgb(j,:))
%         
%         if i<5
%         plot(tExp,yEx(:,i),'.','color',rgb(j,:),'MarkerSize',20)
%         end
%         
%         ylabel(y_axis(i))
%         xlabel('Time (h)')
%         xlim([t(1) t(end)])
%         title(graph_titles{i})
%         set(gca,'Layer','top','XGrid','on','YGrid','on');
%     end
% end
end

