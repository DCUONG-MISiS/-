clc; clear all; close all;

Miu=0.15;
X=[0.5:0.1:4];
for i=1:length(X)
    Del_Del0_MN(i) = fDel_Del0_MN(Miu, X(i));
    Del_Del0_QT(i) = fDel_Del0_QT(Miu, X(i));
    dDel_Del0(i)=Del_Del0_MN(i)-Del_Del0_QT(i);
end

Eps=[0.1:0.01:1.2];
for i=1:length(Eps)
    SigItb_MN(i) = fSigItb_MN(Miu, Eps(i));
    Eta_MN(i) = fEta_MN(Miu, Eps(i));
    F_kNmm_MN(i) = fF_kNmm_MN(Miu, Eps(i));
    SigItb_QT(i) = fSigItb_QT(Miu, Eps(i));
    Eta_QT(i) = fEta_QT(Miu, Eps(i));
    F_kNmm_QT(i) = fF_kNmm_QT(Miu, Eps(i));
    
    dSigItb(i)=SigItb_MN(i)-SigItb_QT(i);
    dEta(i)=Eta_MN(i)-Eta_QT(i);
    dF_kNmm(i)=F_kNmm_MN(i)-F_kNmm_QT(i);
end

figure(1);
subplot(2,1,1);
plot (X,Del_Del0_MN,'k','LineWidth',1.5);
hold on;
plot (X,Del_Del0_QT,'b','LineWidth',1.5);
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
xlabel('$1000\frac{\Delta_{0}}{B_{0}}$ [-]','Interpreter','latex','FontSize',16);
ylabel('$\frac{\Delta}{\Delta_{0}}$ [-]','Interpreter','latex','FontSize',16);
legend({'‘орма крыши','‘орма €блока'},'FontSize',12);
grid on;
grid minor;

subplot(2,1,2);
plot (Eps,Eta_MN,'k','LineWidth',1.5);
hold on;
plot (Eps,Eta_QT,'b','LineWidth',1.5);
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
xlabel('\epsilon [%]','FontSize',16);
ylabel('\xi [-]','FontSize',16);
legend({'‘орма крыши','‘орма €блока'},'FontSize',12);
grid on;
grid minor;

% subplot(2,2,2);
% plot (Eps,SigItb_MN,'k');
% hold on;
% plot (Eps,SigItb_QT,'b');
% ax = gca;
% ax.XAxis.FontSize = 12;
% ax.YAxis.FontSize = 12;
% xlabel('\epsilon [%]','FontSize',16);
% ylabel('\sigma_{itb} [ћѕа]','FontSize',16);
% legend({'MN','QT'},'FontSize',12);
% grid on;
% grid minor;
% 
% subplot(2,2,4);
% plot (Eps,F_kNmm_MN,'k');
% hold on;
% plot (Eps,F_kNmm_QT,'b');
% ax = gca;
% ax.XAxis.FontSize = 12;
% ax.YAxis.FontSize = 12;
% xlabel('\epsilon [%]','FontSize',16);
% ylabel('F [кЌ/мм]','FontSize',16);
% legend({'MN','QT'},'FontSize',12);
% grid on;
% grid minor;
% 
% figure(2);
% subplot(2,2,1);
% plot (X,dDel_Del0,'k','LineWidth',1.5);
% ax = gca;
% ax.XAxis.FontSize = 12;
% ax.YAxis.FontSize = 12;
% xlabel('$1000\frac{\Delta_{0}}{B_{0}}$ [-]','Interpreter','latex','FontSize',16);
% ylabel('d\Delta [-]','FontSize',16);
% grid on;
% grid minor;
% 
% subplot(2,2,2);
% plot (Eps,dSigItb,'b','LineWidth',1.5);
% ax = gca;
% ax.XAxis.FontSize = 12;
% ax.YAxis.FontSize = 12;
% xlabel('\epsilon [%]','FontSize',16);
% ylabel('d\sigma_{itb} [ћѕа]','FontSize',16);
% grid on;
% grid minor;
% 
% subplot(2,2,3);
% plot (Eps,dEta,'b','LineWidth',1.5);
% ax = gca;
% ax.XAxis.FontSize = 12;
% ax.YAxis.FontSize = 12;
% xlabel('\epsilon [%]','FontSize',16);
% ylabel('d\xi [-]','FontSize',16);
% grid on;
% grid minor;
% 
% subplot(2,2,4);
% plot (Eps,dF_kNmm,'b','LineWidth',1.5);
% ax = gca;
% ax.XAxis.FontSize = 12;
% ax.YAxis.FontSize = 12;
% xlabel('\epsilon [%]','FontSize',16);
% ylabel('dF [кЌ/мм]','FontSize',16);
% grid on;
% grid minor;

function Del_Del0_MN = fDel_Del0_MN(Miu, X)
    a1=0.68406;
    a2=0.79233;
    a3=0.84564;
    a4=-0.83249;
%     X=1000*Del0/B0;
    Del_Del0_MN=(a1*Miu+a2)*a3*X^a4;
end

function SigItb_MN = fSigItb_MN(Miu, Eps)
    a1=-0.63122;
    a2=42.65182;
    a3=11.87504;
    a4=0.05689;
%     Eps = fEps(Miu, Del0, B0);
    SigItb_MN=(a1*Miu+a2)*a3*Eps^a4;
end

function Eta_MN = fEta_MN(Miu, Eps)
    a1=0.06921;
    a2=0.0169;
    a3=0.07783;
    a4=-2.48165;
%     Eps = fEps(Miu, Del0, B0);
    Eta_MN=(a1*Miu+a2)*a3*Eps^a4;
end

function F_kNmm_MN = fF_kNmm_MN(Miu, Eps)
    a1=0.30406;
    a2=3.03938;
    a3=2.13434;
    a4=0.04783;
%     Eps = fEps(Miu, Del0, B0);
    F_kNmm_MN=(a1*Miu+a2)*a3*Eps^a4;
end

function Del_Del0_QT = fDel_Del0_QT(Miu, X)
    a1=0.57145;
    a2=0.7906;
    a3=0.83623;
    a4=-1.14359;
%     X=1000*Del0/B0;
    Del_Del0_QT=(a1*Miu+a2)*a3*X^a4;
end

function SigItb_QT = fSigItb_QT(Miu, Eps)
    a1=-0.51298;
    a2=37.62442;
    a3=13.4568;
    a4=0.05402;
%     Eps = fEps(Miu, Del0, B0);
    SigItb_QT=(a1*Miu+a2)*a3*Eps^a4;
end

function Eta_QT = fEta_QT(Miu, Eps)
    a1=0.0158;
    a2=0.01882;
    a3=0.06131;
    a4=-1.3451;
%     Eps = fEps(Miu, Del0, B0);
    Eta_QT=(a1*Miu+a2)*a3*Eps^a4;
end

function F_kNmm_QT = fF_kNmm_QT(Miu, Eps)
    a1=0.25908;
    a2=2.84914;
    a3=2.27735;
    a4=0.04884;
%     Eps = fEps(Miu, Del0, B0);
    F_kNmm_QT=(a1*Miu+a2)*a3*Eps^a4;
end