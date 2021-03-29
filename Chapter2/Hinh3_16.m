clc; clear all; close all;
EpsP=[0:0.001:0.25];
for i=1:length(EpsP)
    SigmaS(i) = fSigmaS(EpsP(i));
end

plot(EpsP,SigmaS,'LineWidth',2);
xticks([0:0.05:0.30]);
xlim([0 0.30]);
ylim([0 160]);
xlabel('\epsilon_{p} [-]','FontWeight','bold');
ylabel('\sigma_{s} [MPa]','FontWeight','bold');
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
grid on;
grid minor;

function SigmaS = fSigmaS(EpsP)
    if EpsP<=0.002
        SigmaS=60;
    else
        SigmaS=194.8*EpsP^0.1895;
    end
end
