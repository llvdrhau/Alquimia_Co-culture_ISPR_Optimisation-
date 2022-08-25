function [] = f_plot_montecarlo(Matrix1,Matrix2,tExp,yExp)
%% This script plots the solution space determined in the Monte Carlo method
% Alberte Regueira
% Department of Chemical Engineering, USC
% November 2019

t=linspace(0,tExp(end),50);
t = union(t,tExp);


matrix1_glu = Matrix1.glu;
matrix1_lac = Matrix1.lac;
matrix1_ye = Matrix1.ye;
matrix1_bm = Matrix1.bm_su;
matrix1_tan = Matrix1.tan;

matrix2_glu = Matrix2.glu;
matrix2_lac = Matrix2.lac;
matrix2_ye = Matrix2.ye;
matrix2_bm = Matrix2.bm_su;
matrix2_tan = Matrix2.tan;


compounds = {'glu','lac','bm','ye'}; %,'tan'};
y_axis = {'Glucose (gCOD/L)','Lactate(gCOD/L)','Biomass (gCOD/L)','Yeast extract (gCOD/L)','mol N/L'};
% graph_titles = {'Glucose','Lactate','Biomass','Yeast extract','TAN'};
% graph characteristics 
%mkSize = 8;
set(gcf,'color','w', 'Position',  [300, 300, 1000, 700]); %background color
%set(gca,'fontsize',14)
set(0,'DefaultLineLineWidth',1,'DefaultAxesFontSize',20,...
    'DefaultAxesXGrid','off','DefaultAxesYGrid','off',...
    'DefaultLineMarkerSize',13.5,...
    'DefaultAxesFontName', 'Calibri', ...  %'Century Schoolbook',...
    'DefaultAxesXcolor', 0.25*[1,1,1],'DefaultAxesYcolor', 0.25*[1,1,1]);

rgb = [186 138 18; 46 76 58; 132 132 132]/255;
rgb_area = [255 220 122; 157 212 179; 179 179 179]/255;

%% letters of subplots
nIDs = 4;
alphabet = ('a':'z').';
chars = num2cell(alphabet(nIDs+1:nIDs+3));
chars = chars.';
charlbl = strcat('(',chars,')'); % {'(a)','(b)','(c)','(d)'}

%% plots 
j =2;
    for i=1:length(compounds)-1
        
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
        
        %plot data
        if i<4
        plot(tExp,yExp(:,i),'^','color',rgb(1,:),'MarkerSize',8,'MarkerFaceColor',rgb(1,:),'DisplayName','dataset 1.1')
        hold on 
        plot(tExp,yExp(:,3+i),'*','color',rgb(2,:),'MarkerSize',10,'DisplayName','dataset 1.2')
        %legend('Monte Carlo confidance interval','','Simulation','Dataset 1.1','Dataset 1.2')
        %legend('','','','Dataset 1.1','Dataset 1.2')
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
    %end
end

