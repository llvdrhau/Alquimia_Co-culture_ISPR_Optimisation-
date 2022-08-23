function [Qgas] = plot_reactor_results(tOut,y, parameters,colour,legendText)
%% Plot results of each reactor simulation

% Alberte Regueira. University of Santiago de Compostela. Spain
% November 2019. Please contact alberte.regueira@gmail.com if you
% intend to use this code.

if nargin<4
    colour='b'; %By default blue
end

%Set up default plot caracteristics
set(0,'DefaultLineLineWidth',1,'DefaultAxesFontSize',11,...
    'DefaultAxesXGrid','off','DefaultAxesYGrid','off',...
    'DefaultLineMarkerSize',4.5,...
    'DefaultAxesFontName','Century Schoolbook',...
    'DefaultAxesXcolor', 0.25*[1,1,1],'DefaultAxesYcolor', 0.25*[1,1,1]);

compounds=parameters.compounds;
indexCompounds=parameters.indexCompounds;
h =  findobj('type','figure');
nPreviousFig = length(h);
nPlots=length(indexCompounds);  %As many plots as states
nFigures=ceil(nPlots/4);        %Plots grouped in 4 subplots figures
iCompound=1;
for i=nPreviousFig+1:nFigures+nPreviousFig
    figure(i)
    for j=1:4
        if iCompound>nPlots     
            break
        end
        subplot(2,2,j)
        yVector=y(:,iCompound); %Select state
        plot(tOut,yVector,colour)
        hold on
        xlabel('Time (d)')
        if indexCompounds(iCompound)==1 %In case it is the total mass
            yLabelText='Total mass (kg)';
        else
            yLabelText=horzcat(compounds{iCompound},' COD kg/m^3');
        end
        ylabel(yLabelText)
        iCompound=iCompound+1;
    end
    if nargin>4
        legend(legendText)
    end
end

figure(nFigures+nPreviousFig+1)
for i=1:length(tOut)
    states=y(i,:)';
    pH(i)=f_pH(states,parameters);
end
plot(tOut,pH)
xlabel('Time (h)')
ylabel('pH')

%% Plot inhibition
for i=1:length(pH)
    pos_Sh2 = find(strcmp(parameters.compoundAbb,'Sh2'));
    Sh2=y(i,pos_Sh2);
    pos_Stic= find(strcmp(parameters.compoundAbb,'Stic'));
    Stic=y(i,pos_Stic);
    inh(i,:)=f_inhibition(pH(i),Sh2,Stic,parameters);
end

inh_names={'Inh. pH (acidogenesis and acetogenesis)';'Inh. pH (acetoclastic meth)';'Inh. pH (hydrogenotrophic meth)';'Inh. H2 (LCFA uptake)';'Inh. H2 (val and but uptake)';'Inh. H2 (prop uptake)';'Monod for In. carbon'};

nPlots=(size(inh,2)/3);
%resto=(nPlots-floor(nPlots))*3; %todo acabar de xeralizar!
iInh=1;
for i=1:nPlots
    figure(nFigures+nPreviousFig+1+i)
    for j=1:3
        subplot(1,3,j)
        plot(tOut,inh(:,iInh),colour)
        hold on
        xlabel('Time (h)')
        ylabel(inh_names{iInh})
        iInh=iInh+1;
    end
    
end
figure
plot(tOut,inh(:,iInh),colour)
xlabel('Time (h)')
ylabel(inh_names{iInh})


for i=1:length(tOut)
    states=y(i,:)';
    Qgas(i)=f_gasPhase(states,parameters);
end

QgasEE=Qgas(end);
fprintf('Gas flow rate is (m3/d): %f \n',QgasEE); 


end










