function Batdau()
    close all; clc; clear all;

    [B0_m Del0_m Miu_m Del_m Eps_m SigRsurf_m SigIsurf_m SigItb_m mSig_m Eta_m F_m] = fLoadDataFromFile('QT.txt');

    X_m=1000*Del0_m./B0_m;
    Del_Del0_m=Del_m./Del0_m;

    N=50;
    Xmin=min(X_m);
    Xmax=max(X_m);
    X=[Xmin:(Xmax-Xmin)/N:Xmax];
    Miu=[0.05:0.01:0.25];
    
    %% 1.////
    for i=1:N+1
        for j=1:length(Miu)
            Del_Del0(i,j) = fDel_Del0(Miu(j), X(i));
        end
    end

    figure(1);
    scatter3(Miu_m,X_m,Del_Del0_m,60,'^','filled');
    hold on;
    surf(Miu,X,Del_Del0,'EdgeColor','none');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.XAxis.Color = 'b';
    ax.YAxis.FontSize = 12;
    ax.YAxis.Color = 'k';
    ax.ZAxis.FontSize = 12;
    ax.ZAxis.Color = 'r';
    xlabel('\mu [-]','FontSize',16);
    ylabel('$1000\frac{\Delta_{0}}{B_{0}}$ [-]','Interpreter','latex','FontSize',16);
    zlabel('$\frac{\Delta}{\Delta_{0}}$ [-]','Interpreter','latex','FontSize',16);
    legend({'������ �������������','�������������'},'FontSize',12);
    grid on;
    grid minor; 

    %% 2////
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
    xlabel('B_{0} [��]','FontSize',16);
    ylabel('\epsilon [%]','FontSize',16);
    legend({'������ �������������','�������������'},'FontSize',12);
    grid on;
    grid minor;

    %% 3///////
    Del0_SigRes_onSurf_m=Del0_m.*SigRsurf_m;
    for i=1:N+1
        for j=1:length(Miu)
            Del0_SigRes_onSurf(i,j) = fDel0_SigRes_onSurf(Miu(j), X(i));
        end
    end

    figure(3);
    scatter3(Miu_m,X_m,Del0_SigRes_onSurf_m,60,'^','filled');
    hold on;
    surf(Miu,X,Del0_SigRes_onSurf,'EdgeColor','none');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.XAxis.Color = 'b';
    ax.YAxis.FontSize = 12;
    ax.YAxis.Color = 'k';
    ax.ZAxis.FontSize = 12;
    ax.ZAxis.Color = 'r';
    xlabel('\mu [-]','FontSize',16);
    ylabel('$1000\frac{\Delta_{0}}{B_{0}}$ [-]','Interpreter','latex','FontSize',16);
    zlabel('\Delta_{0}.\sigma_{rsurf} [���.��]','FontSize',16);
    legend({'������ �������������','�������������'},'FontSize',12);
    grid on;
    grid minor;

    %% 4////
    Del0_SigI_onSurf_m=Del0_m.*SigIsurf_m;

    for i=1:N+1
        for j=1:length(Miu)
            Del0_SigI_onSurf(i,j) = fDel0_SigI_onSurf(Miu(j), X(i));
        end
    end

    figure(4);
    scatter3(Miu_m,X_m,Del0_SigI_onSurf_m,60,'^','filled');
    hold on;
    surf(Miu,X,Del0_SigI_onSurf,'EdgeColor','none');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.XAxis.Color = 'b';
    ax.YAxis.FontSize = 12;
    ax.YAxis.Color = 'k';
    ax.ZAxis.FontSize = 12;
    ax.ZAxis.Color = 'r';
    xlabel('\mu [-]','FontSize',16);
    ylabel('$1000\frac{\Delta_{0}}{B_{0}}$ [-]','Interpreter','latex','FontSize',16);
    zlabel('\Delta_{0}.\sigma_{isurf} [���.��]','FontSize',16);
    legend({'������ �������������','�������������'},'FontSize',12);
    grid on;
    grid minor;

    %% 5////
    for i=1:N+1
        for j=1:length(Miu)
            SigItb(i,j) = fSigItb(Miu(j), Eps(i));
        end
    end

    figure(5);
    scatter3(Miu_m,Eps_m,SigItb_m,60,'^','filled');
    hold on;
    surf(Miu,Eps,SigItb,'EdgeColor','none');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.XAxis.Color = 'b';
    ax.YAxis.FontSize = 12;
    ax.YAxis.Color = 'k';
    ax.ZAxis.FontSize = 12;
    ax.ZAxis.Color = 'r';
    xlabel('\mu [-]','FontSize',16);
    ylabel('\epsilon [%]','FontSize',16);
    zlabel('\sigma_{itb} [���]','FontSize',16);
    legend({'������ �������������','�������������'},'FontSize',12);
    grid on;
    grid minor;

    %% 6///////
    Miu_mSig_m=Miu_m.*mSig_m;
    for i=1:N+1
        for j=1:length(Miu)
            Miu_mSig(i,j) = fMiu_mSig(Miu(j), X(i));
        end
    end

    figure(6);
    scatter3(Miu_m,X_m,Miu_mSig_m,60,'^','filled');
    hold on;
    surf(Miu,X,Miu_mSig,'EdgeColor','none');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.XAxis.Color = 'b';
    ax.YAxis.FontSize = 12;
    ax.YAxis.Color = 'k';
    ax.ZAxis.FontSize = 12;
    ax.ZAxis.Color = 'r';
    xlabel('\mu [-]','FontSize',16);
    ylabel('$1000\frac{\Delta_{0}}{B_{0}}$ [-]','Interpreter','latex','FontSize',16);
    zlabel('\mu.m_{\sigma} [-]','FontSize',16);
    legend({'������ �������������','�������������'},'FontSize',12);
    grid on;
    grid minor;

    %% 7////
    for i=1:N+1
        for j=1:length(Miu)
            Eta(i,j) = fEta(Miu(j), Eps(i));
        end
    end

    figure(7);
    scatter3(Miu_m,Eps_m,Eta_m,60,'^','filled');
    hold on;
    surf(Miu,Eps,Eta,'EdgeColor','none');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.XAxis.Color = 'b';
    ax.YAxis.FontSize = 12;
    ax.YAxis.Color = 'k';
    ax.ZAxis.FontSize = 12;
    ax.ZAxis.Color = 'r';
    xlabel('\mu [-]','FontSize',16);
    ylabel('\epsilon [%]','FontSize',16);
    zlabel('\xi [-]','FontSize',16);
    legend({'������ �������������','�������������'},'FontSize',12);
    grid on;
    grid minor;

    %% 8////
    F_kNmm_m=F_m/45;
    for i=1:N+1
        for j=1:length(Miu)
            F_kNmm(i,j) = fF_kNmm(Miu(j), Eps(i));
        end
    end

    figure(8);
    scatter3(Miu_m,Eps_m,F_kNmm_m*2.849,60,'^','filled');
    hold on;
    surf(Miu,Eps,F_kNmm*2.849,'EdgeColor','none');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.XAxis.Color = 'b';
    ax.YAxis.FontSize = 12;
    ax.YAxis.Color = 'k';
    ax.ZAxis.FontSize = 12;
    ax.ZAxis.Color = 'r';
    xlabel('\mu [-]','FontSize',16);
    ylabel('\epsilon [%]','FontSize',16);
    zlabel('p [���]','FontSize',16);
    legend({'������ �������������','�������������'},'FontSize',12);
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