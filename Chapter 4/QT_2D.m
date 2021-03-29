function Batdau()
    close all; clc; clear all;

    [B0_m1 Del0_m1 Miu_m1 Del_m1 Eps_m1 SigRsurf_m1 SigIsurf_m1 SigItb_m1 mSig_m1 Eta_m1 F_m1] = fLoadDataFromFile('QT_005.txt');
    [B0_m2 Del0_m2 Miu_m2 Del_m2 Eps_m2 SigRsurf_m2 SigIsurf_m2 SigItb_m2 mSig_m2 Eta_m2 F_m2] = fLoadDataFromFile('QT_015.txt');
    [B0_m3 Del0_m3 Miu_m3 Del_m3 Eps_m3 SigRsurf_m3 SigIsurf_m3 SigItb_m3 mSig_m3 Eta_m3 F_m3] = fLoadDataFromFile('QT_025.txt');
    [B0_m Del0_m Miu_m Del_m Eps_m SigRsurf_m SigIsurf_m SigItb_m mSig_m Eta_m F_m] = fLoadDataFromFile('QT.txt');

    X_m=1000*Del0_m1./B0_m1;
    Del_Del0_m1=Del_m1./Del0_m1;
    Del_Del0_m2=Del_m2./Del0_m2;
    Del_Del0_m3=Del_m3./Del0_m3;

    N=50;
    Xmin=min(X_m);
    Xmax=max(X_m);
    X=[Xmin:(Xmax-Xmin)/N:Xmax];
    Miu1=0.05;
    Miu2=0.15;
    Miu3=0.25;
    
    %% 1.//////
    for i=1:N+1
        Del_Del01(i) = fDel_Del0(Miu1, X(i));
        Del_Del02(i) = fDel_Del0(Miu2, X(i));
        Del_Del03(i) = fDel_Del0(Miu3, X(i));
    end

    figure(1);
    scatter(X_m,Del_Del0_m1,60,'^','filled');
    hold on;
    plot(X,Del_Del01,'LineWidth',1.5);
    hold on;
    scatter(X_m,Del_Del0_m2,60,'s','filled');
    hold on;
    plot(X,Del_Del02,'LineWidth',1.5);
    hold on;
    scatter(X_m,Del_Del0_m3,60,'d','filled');
    hold on;
    plot(X,Del_Del03,'LineWidth',1.5);
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    xlabel('$1000\frac{\Delta_{0}}{B_{0}}$ [-]','Interpreter','latex','FontSize',16);
    ylabel('$\frac{\Delta}{\Delta_{0}}$ [-]','Interpreter','latex','FontSize',16);
    legend({'Данные при \mu=0,05','Аппроксимация при \mu=0,05','Данные при \mu=0,15','Аппроксимация при \mu=0,15','Данные при \mu=0,25','Аппроксимация при \mu=0,25'},'FontSize',12);
    grid on;
    grid minor; 
    
    %% 2./////
    B0min=min(B0_m);
    B0max=max(B0_m);
    B0=[B0min:(B0max-B0min)/N:B0max];
    Eps = fEps(B0);
    figure(2);
    scatter(B0_m,Eps_m,60,'^');
    hold on;
    plot(B0,Eps,'LineWidth',1.5);
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    xlabel('B_{0} [мм]','FontSize',16);
    ylabel('\epsilon [%]','FontSize',16);
    legend({'Данные моделирования','Аппроксимация'},'FontSize',12);
    grid on;
    grid minor;

    %% 3///////
    Del0_SigRes_onSurf_m1=Del0_m1.*SigRsurf_m1;
    Del0_SigRes_onSurf_m2=Del0_m2.*SigRsurf_m2;
    Del0_SigRes_onSurf_m3=Del0_m3.*SigRsurf_m3;
    for i=1:N+1
        Del0_SigRes_onSurf1(i) = fDel0_SigRes_onSurf(Miu1, X(i));
        Del0_SigRes_onSurf2(i) = fDel0_SigRes_onSurf(Miu2, X(i));
        Del0_SigRes_onSurf3(i) = fDel0_SigRes_onSurf(Miu3, X(i));
    end

    figure(3);
    scatter(X_m,Del0_SigRes_onSurf_m1,60,'^','filled');
    hold on;
    plot(X,Del0_SigRes_onSurf1,'LineWidth',1.5);
    hold on;
    scatter(X_m,Del0_SigRes_onSurf_m2,60,'s','filled');
    hold on;
    plot(X,Del0_SigRes_onSurf2,'LineWidth',1.5);
    hold on;
    scatter(X_m,Del0_SigRes_onSurf_m3,60,'d','filled');
    hold on;
    plot(X,Del0_SigRes_onSurf3,'LineWidth',1.5);
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    xlabel('$1000\frac{\Delta_{0}}{B_{0}}$ [-]','Interpreter','latex','FontSize',16);
    ylabel('\Delta_{0}.\sigma_{rsurf} [МПа.мм]','FontSize',16);
    legend({'Данные при \mu=0,05','Аппроксимация при \mu=0,05','Данные при \mu=0,15','Аппроксимация при \mu=0,15','Данные при \mu=0,25','Аппроксимация при \mu=0,25'},'FontSize',12);
    grid on;
    grid minor;

    %% 4////
    Del0_SigI_onSurf_m1=Del0_m1.*SigIsurf_m1;
    Del0_SigI_onSurf_m2=Del0_m2.*SigIsurf_m2;
    Del0_SigI_onSurf_m3=Del0_m3.*SigIsurf_m3;
    for i=1:N+1
        Del0_SigI_onSurf1(i) = fDel0_SigI_onSurf(Miu1, X(i));
        Del0_SigI_onSurf2(i) = fDel0_SigI_onSurf(Miu2, X(i));
        Del0_SigI_onSurf3(i) = fDel0_SigI_onSurf(Miu3, X(i));
    end

    figure(4);
    scatter(X_m,Del0_SigI_onSurf_m1,60,'^','filled');
    hold on;
    plot(X,Del0_SigI_onSurf1,'LineWidth',1.5);
    hold on;
    scatter(X_m,Del0_SigI_onSurf_m2,60,'s','filled');
    hold on;
    plot(X,Del0_SigI_onSurf2,'LineWidth',1.5);
    hold on;
    scatter(X_m,Del0_SigI_onSurf_m3,60,'d','filled');
    hold on;
    plot(X,Del0_SigI_onSurf3,'LineWidth',1.5);
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    xlabel('$1000\frac{\Delta_{0}}{B_{0}}$ [-]','Interpreter','latex','FontSize',16);
    ylabel('\Delta_{0}.\sigma_{isurf} [МПа.мм]','FontSize',16);
    legend({'Данные при \mu=0,05','Аппроксимация при \mu=0,05','Данные при \mu=0,15','Аппроксимация при \mu=0,15','Данные при \mu=0,25','Аппроксимация при \mu=0,25'},'FontSize',12);
    grid on;
    grid minor;

    %% 5////
    for i=1:N+1
        SigItb1(i) = fSigItb(Miu1, Eps(i));
        SigItb2(i) = fSigItb(Miu2, Eps(i));
        SigItb3(i) = fSigItb(Miu3, Eps(i));
    end

    figure(5);
    scatter(Eps_m1,SigItb_m1,60,'^','filled');
    hold on;
    plot(Eps,SigItb1,'LineWidth',1.5);
    hold on;
    scatter(Eps_m2,SigItb_m2,60,'s','filled');
    hold on;
    plot(Eps,SigItb2,'LineWidth',1.5);
    hold on;
    scatter(Eps_m3,SigItb_m3,60,'d','filled');
    hold on;
    plot(Eps,SigItb3,'LineWidth',1.5);
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    xlabel('\epsilon [%]','FontSize',16);
    ylabel('\sigma_{itb} [МПа]','FontSize',16);
    legend({'Данные при \mu=0,05','Аппроксимация при \mu=0,05','Данные при \mu=0,15','Аппроксимация при \mu=0,15','Данные при \mu=0,25','Аппроксимация при \mu=0,25'},'FontSize',12);
    grid on;
    grid minor;

    %% 6///////
    Miu_mSig_m1=Miu_m1.*mSig_m1;
    Miu_mSig_m2=Miu_m2.*mSig_m2;
    Miu_mSig_m3=Miu_m3.*mSig_m3;
    for i=1:N+1
        Miu_mSig1(i) = fMiu_mSig(Miu1, X(i));
        Miu_mSig2(i) = fMiu_mSig(Miu2, X(i));
        Miu_mSig3(i) = fMiu_mSig(Miu3, X(i));
    end

    figure(6);
    scatter(X_m,Miu_mSig_m1,60,'^','filled');
    hold on;
    plot(X,Miu_mSig1,'LineWidth',1.5);
    hold on;
    scatter(X_m,Miu_mSig_m2,60,'s','filled');
    hold on;
    plot(X,Miu_mSig2,'LineWidth',1.5);
    hold on;
    scatter(X_m,Miu_mSig_m3,60,'d','filled');
    hold on;
    plot(X,Miu_mSig3,'LineWidth',1.5);
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    xlabel('$1000\frac{\Delta_{0}}{B_{0}}$ [-]','Interpreter','latex','FontSize',16);
    ylabel('\mu.m_{\sigma} [-]','FontSize',16);
    legend({'Данные при \mu=0,05','Аппроксимация при \mu=0,05','Данные при \mu=0,15','Аппроксимация при \mu=0,15','Данные при \mu=0,25','Аппроксимация при \mu=0,25'},'FontSize',12);
    grid on;
    grid minor;

    %% 7////
    for i=1:N+1
        Eta1(i) = fEta(Miu1, Eps(i));
        Eta2(i) = fEta(Miu2, Eps(i));
        Eta3(i) = fEta(Miu3, Eps(i));
    end

    figure(7);
    scatter(Eps_m1,Eta_m1,60,'^','filled');
    hold on;
    plot(Eps,Eta1,'LineWidth',1.5);
    hold on;
    scatter(Eps_m2,Eta_m2,60,'s','filled');
    hold on;
    plot(Eps,Eta2,'LineWidth',1.5);
    hold on;
    scatter(Eps_m3,Eta_m3,60,'d','filled');
    hold on;
    plot(Eps,Eta3,'LineWidth',1.5);
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    xlabel('\epsilon [%]','FontSize',16);
    ylabel('\xi [-]','FontSize',16);
    legend({'Данные при \mu=0,05','Аппроксимация при \mu=0,05','Данные при \mu=0,15','Аппроксимация при \mu=0,15','Данные при \mu=0,25','Аппроксимация при \mu=0,25'},'FontSize',12);
    grid on;
    grid minor;

    %% 8////
    F_kNmm_m1=F_m1/45;
    F_kNmm_m2=F_m2/45;
    F_kNmm_m3=F_m3/45;
    for i=1:N+1
        F_kNmm1(i) = fF_kNmm(Miu1, Eps(i));
        F_kNmm2(i) = fF_kNmm(Miu2, Eps(i));
        F_kNmm3(i) = fF_kNmm(Miu3, Eps(i));
    end

    figure(8);
    scatter(Eps_m1,F_kNmm_m1*2.849,60,'^','filled');
    hold on;
    plot(Eps,F_kNmm1*2.849,'LineWidth',1.5);
    hold on;
    scatter(Eps_m2,F_kNmm_m2*2.849,60,'s','filled');
    hold on;
    plot(Eps,F_kNmm2*2.849,'LineWidth',1.5);
    hold on;
    scatter(Eps_m3,F_kNmm_m3*2.849,60,'d','filled');
    hold on;
    plot(Eps,F_kNmm3*2.849,'LineWidth',1.5);
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    xlabel('\epsilon [%]','FontSize',16);
    ylabel('p [МПа]','FontSize',16);
    legend({'Данные при \mu=0,05','Аппроксимация при \mu=0,05','Данные при \mu=0,15','Аппроксимация при \mu=0,15','Данные при \mu=0,25','Аппроксимация при \mu=0,25'},'FontSize',12);
    grid on;
    grid minor;
end

function [B0_m Del0_m Miu_m Del_m Eps_m SigRsurf_m SigIsurf_m SigItb_m mSig_m Eta_m F_m] = fLoadDataFromFile(fileName)
    %% LAY SO LIEU TU FILE
    fileID = fopen(fileName, 'r');
    A=fscanf(fileID,'%f %f',[11 Inf]);
    fclose(fileID);
    B=A';

    B0_m=B(:,1);
    Del0_m=B(:,2);
    Miu_m=B(:,3);
    Del_m=B(:,4);
    Eps_m=B(:,5);
    SigRsurf_m=B(:,6);
    SigIsurf_m=B(:,7);
    SigItb_m=B(:,8);
    mSig_m=B(:,9);
    Eta_m=B(:,10);
    F_m=B(:,11);
end

%% 1. Do oval (mm)
function Del_Del0 = fDel_Del0(Miu, X)
    a1=0.57145;
    a2=0.7906;
    a3=0.83623;
    a4=-1.14359;
%     X=1000*Del0/B0;
    Del_Del0=(a1*Miu+a2)*a3*X^a4;
end

%% 2. Muc do bien dang (%)
function Eps = fEps(B0)
    a1=-21.8193;
    a2=95.87457;
    Eps=a1*B0/1000+a2;
end

%% 3. Ung suat du lon nhat tren be mat ngoai (MPa)
function Del0_SigRes_onSurf = fDel0_SigRes_onSurf(Miu, X) 
    a1=0.14068;
    a2=1.63283;
    a3=2.2076;
    a4=0.45385;
%     X=1000*Del0/B0;
    Del0_SigRes_onSurf=100*(a1*Miu+a2)*(a3*X+a4);
end

%% 4. Cuong do ung suat lon nhat tren be mat ngoai (MPa)
function Del0_SigI_onSurf = fDel0_SigI_onSurf(Miu, X)
    a1=-0.0795;
    a2=4.22179;
    a3=5.79337;
    a4=-0.33984;
%     X=1000*Del0/B0;
    Del0_SigI_onSurf=100*(a1*Miu+a2)*(a3*X+a4);
end

%% 5. Cuong do ung suat trung binh cua toan vat the (MPa)
function SigItb = fSigItb(Miu, Eps)
    a1=-0.51298;
    a2=37.62442;
    a3=13.4568;
    a4=0.05402;
%     Eps = fEps(Miu, Del0, B0);
    SigItb=(a1*Miu+a2)*a3*Eps^a4;
end

%% 6. He so tap trung ung suat (-)
function Miu_mSig = fMiu_mSig(Miu, X)
    a1=2.38089;
    a2=-0.000875969;
    a3=0.0092;
    a4=0.43565;
%     X=1000*Del0/B0;
    Miu_mSig=(a1*Miu+a2)*(a3*X+a4);
end

%% 7. He so dan hoi du (-)
function Eta = fEta(Miu, Eps)
    a1=0.0158;
    a2=0.01882;
    a3=0.06131;
    a4=-1.3451;
%     Eps = fEps(Miu, Del0, B0);
    Eta=(a1*Miu+a2)*a3*Eps^a4;
end

%% 8. Luc lon nhat cua 1 Segment vao ong cho 1mm chieu dai gian no (kN/mm)
function F_kNmm = fF_kNmm(Miu, Eps)
    a1=0.25908;
    a2=2.84914;
    a3=2.27735;
    a4=0.04884;
%     Eps = fEps(Miu, Del0, B0);
    F_kNmm=(a1*Miu+a2)*a3*Eps^a4;
end