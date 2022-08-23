function inh=f_inhibition(pH,Sh2,Stic,parameters)
%% Determine inhibition values at each time point

% Alberte Regueira. University of Santiago de Compostela. Spain
% November 2019. Please contact alberte.regueira@gmail.com if you
% intend to use this code.

parValues = parameters.parValues;
parAbb = parameters.parAbb;

pHUL_aa = parValues(strcmp(parAbb,'pHUL_aa'));
pHLL_aa = parValues(strcmp(parAbb,'pHLL_aa'));
pHUL_h2 = parValues(strcmp(parAbb,'pHUL_h2'));
pHLL_h2 = parValues(strcmp(parAbb,'pHLL_h2'));
pHUL_ac = parValues(strcmp(parAbb,'pHUL_ac'));
pHLL_ac = parValues(strcmp(parAbb,'pHLL_ac'));
KI_h2_pro = parValues(strcmp(parAbb,'KI_h2_pro'));
KI_h2_c4 = parValues(strcmp(parAbb,'KI_h2_c4'));
KI_h2_fa = parValues(strcmp(parAbb,'KI_h2_fa'));

naa = 3.0/(pHUL_aa-pHLL_aa);
nac = 3.0/(pHUL_ac-pHLL_ac);
nh2 = 3.0/(pHUL_h2-pHLL_h2);
KpH_aa=10^-(0.5*(pHUL_aa+pHLL_aa));
KpH_ac=10^-(0.5*(pHUL_ac+pHLL_ac));
KpH_h2=10^-(0.5*(pHUL_h2+pHLL_h2));
sH = 10^(-pH);

IpH_aa = KpH_aa^naa /(sH^naa + KpH_aa^naa);
IpH_ac = KpH_ac^nac /(sH^nac + KpH_ac^nac);
IpH_h2 = KpH_h2^nh2 /(sH^nh2 + KpH_h2^nh2);



inh(1)=IpH_aa;
inh(2)=IpH_ac;
inh(3)=IpH_h2;

Ih2_fa = 1/(1+(Sh2/KI_h2_fa)); 
Ih2_c4 = 1/(1+(Sh2/KI_h2_c4));
Ih2_pro = 1/(1+(Sh2/KI_h2_pro));

inh(4)=Ih2_fa;
inh(5)=Ih2_c4;
inh(6)=Ih2_pro;

Mon_tic=Stic./(3e-4+Stic);
inh(7)=Mon_tic;




end
