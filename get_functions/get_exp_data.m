function [yExp,tExp,data] = get_exp_data(dataName)
% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
%   Retrives experimental data from excel files 
%   
switch dataName
    case 'V3'
        [data1,Names]=xlsread('expData.xlsx','V3');
        tExp=data1(:,strcmp(Names,'time'));
        
        yExp1(:,1)=data1(:,strcmp(Names,'Lac'));
        yExp1(:,2)=data1(:,strcmp(Names,'Pro'));
        yExp1(:,3)=data1(:,strcmp(Names,'Ace'));
        yExp1(:,4)=data1(:,strcmp(Names,'BM'));
        yExp = yExp1;
        
        data = data1;
    case 'V4'
        [data1,Names]=xlsread('expData.xlsx','V4');
        tExp=data1(:,strcmp(Names,'time'));
        
        yExp1(:,1)=data1(:,strcmp(Names,'Lac'));
        yExp1(:,2)=data1(:,strcmp(Names,'Pro'));
        yExp1(:,3)=data1(:,strcmp(Names,'Ace'));
        yExp1(:,4)=data1(:,strcmp(Names,'BM'));
        yExp = yExp1;
        
        data = data1;
    case 'V34'
        [data1,~]=xlsread('expData.xlsx','V3');
        [data2,Names]=xlsread('expData.xlsx','V4');
        
        tExp=data1(:,strcmp(Names,'time'));
        yExp1(:,1)=data1(:,strcmp(Names,'Lac'));
        yExp1(:,2)=data1(:,strcmp(Names,'Pro'));
        yExp1(:,3)=data1(:,strcmp(Names,'Ace'));
        yExp1(:,4)=data1(:,strcmp(Names,'BM'));
        
        % tExp2=data2(:,strcmp(Names,'t'));
        yExp2(:,1)=data2(:,strcmp(Names,'Lac'));
        yExp2(:,2)=data2(:,strcmp(Names,'Pro'));
        yExp2(:,3)=data2(:,strcmp(Names,'Ace'));
        yExp2(:,4)=data2(:,strcmp(Names,'BM'));
        
        yExp = [yExp1 yExp2];
        data = [data1 data2];
    case 'BV5'
        [data1,Names]=xlsread('expData.xlsx','BV5');
        
        tExp=data1(:,strcmp(Names,'time'));
        
        yExp1(:,1)=data1(:,strcmp(Names,'Glu'));
        yExp1(:,2)=data1(:,strcmp(Names,'Lac'));
        yExp1(:,3)=data1(:,strcmp(Names,'Pro'));
        yExp1(:,4)=data1(:,strcmp(Names,'Ace'));
        yExp1(:,5)=data1(:,strcmp(Names,'BM'));
        yExp = yExp1; 
        data = data1;
end    
end

