clc; clear all; close all;
EpsP=[0:0.001:0.25];
for i=1:length(EpsP)
    SigmaS1(i) = fSigmaS1(EpsP(i));
    SigmaS2(i) = fSigmaS2(EpsP(i));
    SigmaS3(i) = fSigmaS3(EpsP(i));
end

figure(1);
subplot(2,2,1);
plot(EpsP,SigmaS1,'LineWidth',2);
xticks([0:0.05:0.30]);
xlim([0 0.30]);
ylim([0 650]);
xlabel('\epsilon_{p} [-]','FontWeight','bold');
ylabel('\sigma_{s} [MPa]','FontWeight','bold');
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
grid on;
grid minor;

subplot(2,2,2);
plot(EpsP,SigmaS2,'LineWidth',2);
xticks([0:0.05:0.30]);
xlim([0 0.30]);
ylim([0 650]);
xlabel('\epsilon_{p} [-]','FontWeight','bold');
ylabel('\sigma_{s} [MPa]','FontWeight','bold');
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
grid on;
grid minor;

subplot(2,2,3);
plot(EpsP,SigmaS3,'LineWidth',2);
xticks([0:0.05:0.30]);
xlim([0 0.30]);
ylim([0 650]);
xlabel('\epsilon_{p} [-]','FontWeight','bold');
ylabel('\sigma_{s} [MPa]','FontWeight','bold');
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
grid on;
grid minor;

function SigmaS1 = fSigmaS1(EpsP)
    SigmaS1=460;
end

function SigmaS2 = fSigmaS2(EpsP)
        SigmaS2=460+659.3*EpsP;
end

function SigmaS3 = fSigmaS3(EpsP)
    if EpsP<=0.002
        SigmaS3=460;
    else
        SigmaS3=644.3*EpsP^0.0542;
    end
end
