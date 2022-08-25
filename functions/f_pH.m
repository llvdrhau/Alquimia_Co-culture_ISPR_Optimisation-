function [pH,residual]=f_pH(states,parameters,pHguess,t)
% Lucas Van der Hauwaert & Alberte Regueira. University of Santiago de 
% Compostela. Spain October 2021.Please contact 
% lucas.vanderhauwaert@usc.es if you intend to use this code.

% Solution of pH equations by combination of the mass balances, equlibrium
% and charge balance into one nonlinear equation. The equation is solved by
% Brent-Dekker's method

% Alberte Regueira. University of Santiago de Compostela. Spain
% November 2019. Please contact alberte.regueira@gmail.com if you
% intend to use this code.

%If the pH is assumed constant, end the function and return a value stored
%in parameters

%If flag pHconstant, return the value given and end the function execution

% In the case of the co_culture project of BioChem pH is assumed cte and
% this script is not used in any significant way!


pHconstant=parameters.pH.pHconstant;
if pHconstant
    if isfield(parameters.pH,'pHtrajectory')    %If pH is not constant during the simulation
        pHvector=parameters.pH.pHtrajectory(:,2);
        tpH=parameters.pH.pHtrajectory(:,1);
        if length(tpH)>1
            pH = interp1(tpH,pHvector,t);
        else
            pH=pHvector;
        end
    else
        pH=parameters.pH.pHValue;
    end
    residual=999;
    return
end

indexCompounds=parameters.indexCompounds;
Mw=parameters.Mw;
parValues = parameters.parValues;
parAbb = parameters.parAbb;
compoundAbb=parameters.compoundAbb;

% Identify compounds
Sva = states(strcmp(compoundAbb,'Sva'));        %Valerate
Sbu = states(strcmp(compoundAbb,'Sbu'));        %Butyrate
Spro = states(strcmp(compoundAbb,'Spro'));      %Propionate
Slac = states(strcmp(compoundAbb,'Slac'));      %Lactate
Sac = states(strcmp(compoundAbb,'Sac'));        %Acetate
Stic = states(strcmp(compoundAbb,'Stic'));      %TIC
SIn = states(strcmp(compoundAbb,'SIn'));        %Inorganic N
Scat = states(strcmp(compoundAbb,'Scat'));      %Cations
San = states(strcmp(compoundAbb,'San'));        %Anions

% Get parameter values
Kac_tic = parValues(strcmp(parAbb,'Kac_tic'));
Kac_ac = parValues(strcmp(parAbb,'Kac_ac'));
Kac_bu = parValues(strcmp(parAbb,'Kac_bu'));
Kac_in = parValues(strcmp(parAbb,'Kac_in'));
Kac_va = parValues(strcmp(parAbb,'Kac_va'));
Kac_pro = parValues(strcmp(parAbb,'Kac_pro'));
Kw = parValues(strcmp(parAbb,'Kw'));
Kac_lac = parValues(strcmp(parAbb,'Kac_lac'));


Z=Scat-San;

%Transform into molar concentrations
Mva = Sva/Mw(strcmp(compoundAbb,'Sva'));
Mbu = Sbu/Mw(strcmp(compoundAbb,'Sbu'));
Mpro = Spro/Mw(strcmp(compoundAbb,'Spro'));
Mac = Sac/Mw(strcmp(compoundAbb,'Sac'));
Mtic = Stic/Mw(strcmp(compoundAbb,'Stic'));
MIn = SIn/Mw(strcmp(compoundAbb,'SIn'));
Mlac= Slac/Mw(strcmp(compoundAbb,'Slac'));

inputPar=[Kac_tic,Kac_ac,Kac_bu,Kac_in,Kac_va,Kac_pro,Kw,Kac_lac,Mva,Mbu,Mpro,Mac,Mtic,MIn,Mlac,Z];

if nargin<3
    pHguess=5;
end

tolX=1e-8;
tolF=1e-8;
deltaMethod=10*tolX;
%% Brent-Dekker's method
%Try to  use a hot start
a=pHguess*1.05;
b=pHguess*0.95;
fa = f_charge_balance(a,inputPar);
fb = f_charge_balance(b,inputPar);
%Check if they are at both sides of the root, if not try farther away from
%the guess
if fa*fb>=0
    a=pHguess*1.5;
    b=pHguess*0.5;
    fa = f_charge_balance(a,inputPar);
    fb = f_charge_balance(b,inputPar);
end
%If they don't bracket the root, set the root to be between pH 0 and 14.
%pH must be between these two values in aqueous solution
if fa*fb>=0
    a=0;
    b=14;
    fa = f_charge_balance(a,inputPar);
    fb = f_charge_balance(b,inputPar);
end

%Set b as the value closest to the root (with the lowest residual)
if abs(fa)<abs(fb)
    old=[a b];
    a=old(2);
    b=old(1);
end
%Initialize c and d
c=a;
d=a;
fa = f_charge_balance(a,inputPar);
fb = f_charge_balance(b,inputPar);
fc = f_charge_balance(c,inputPar);
fs=1;
mflag=1;%mflag can be 1 'set' of 0 'cleared'
nIterations=0;

while abs(fs)>tolF && abs(b-a)>tolX && nIterations<100;
    if fa~=fc && fb~=fc                     % Inverse quadratic interpolation
        A=a*fb*fc/(fa-fb)/(fa-fc);
        B=b*fa*fc/(fb-fa)/(fb-fc);
        C=c*fa*fb/(fc-fa)/(fc-fb);
        s=A+B+C;
    else                                    % Secant method
        s=b-fb*(b-a)/(fb-fa);
    end
    
    %Equivalent code but slower
    %vector1=[(3*a+b)/4 s b];
    %[sortedVector,indexSorted]=sort(vector1);
    % s is not between  (3a+b)/4  and b
    if s>(3*a+b)/4
        if s<b; cond(1)=1; else cond(1)=0;end
    else
        if s<b; cond(1)=0; else cond(1)=1; end
    end
    % Conditions that would mean that the bisection is faster than the other
    % methods
    %cond(1)      s is not between  (3a+b)/4  and b
    cond(2) = (mflag==1&&abs(s-b)>=abs((b-c)/2));   % mflag is set and |s-b| and |b-c|/2) or
    cond(3) = (mflag==0&&abs(s-b)>=abs((c-d)/2));   % mflag is cleared and |s-b| and |c-d|/2) or
    cond(4) = (mflag==1&&abs(b-c)<deltaMethod);     % mflag is set and |b-c| < |delta|) or
    cond(5) = (mflag==0&&abs(c-d)<deltaMethod);     % mflag is cleared and |c-d| < |delta|)
    if any(cond)    %Bisection
        s=(a+b)/2;
        mflag=1;
    else
        mflag=0;
    end
    fs= f_charge_balance(s,inputPar);
    d=c;
    c=b;
    if fa*fs<0 %if f(a) f(s) < 0 then b := s else a := s
        b=s;
        fb = f_charge_balance(b,inputPar);
    else
        a=s;
        fa = f_charge_balance(a,inputPar);
    end
    
    %Set b as the value closest to the root (with the lowest residual)
    if abs(fa)<abs(fb)
        old=[a b fa fb];
        a=old(2);
        b=old(1);
        fa = old(4);
        fb = old(3);
    end
    nIterations=nIterations+1;
end
pH=b;
residual=fb;