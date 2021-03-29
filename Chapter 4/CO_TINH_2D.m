function Batdau()
    close all; clc; clear all;

    [Sig02_m E_1000_m Eps_m LamdaD_m SigRsurf_m] = fLoadDataFromFile('Cotinh.txt');

    X_m=Sig02_m./E_1000_m;
    SigRsurf_E_m=SigRsurf_m./E_1000_m;

    N=50;
    Xmin=min(X_m);
    Xmax=max(X_m);
    X=[Xmin:(Xmax-Xmin)/N:Xmax];
    
    for i=1:N+1
        Eps(i) = fEps(X(i));
        LamdaD(i) = fLamdaD(X(i));
        SigRsurf_E(i) = fSigRsurf_E(X(i));
    end

    %% 1. ///
    figure(1);
    scatter(X_m,Eps_m,60,'^','filled');
    hold on;
    plot(X,Eps,'LineWidth',1.5);
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    xlabel('$1000\frac{\sigma_{0,2}}{E}$ [-]','Interpreter','latex','FontSize',16);
    ylabel('\epsilon [%]','FontSize',16);
    legend({'Данные моделирования','Аппроксимация'},'FontSize',12);
    grid on;
    grid minor; 
    
    %% 2./////
    figure(2);
    scatter(X_m,LamdaD_m,60,'^','filled');
    hold on;
    plot(X,LamdaD,'LineWidth',1.5);
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    xlabel('$1000\frac{\sigma_{0,2}}{E}$ [-]','Interpreter','latex','FontSize',16);
    ylabel('\lambda [%]','FontSize',16);
    legend({'Данные моделирования','Аппроксимация'},'FontSize',12);
    grid on;
    grid minor;

    %% 3./////
    figure(3);
    scatter(X_m,SigRsurf_E_m,60,'^','filled');
    hold on;
    plot(X,SigRsurf_E,'LineWidth',1.5);
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    xlabel('$1000\frac{\sigma_{0,2}}{E}$ [-]','Interpreter','latex','FontSize',16);
    ylabel('$1000\frac{\sigma_{rsurf}}{E}$ [-]','Interpreter','latex','FontSize',16);
    legend({'Данные моделирования','Аппроксимация'},'FontSize',12);
    grid on;
    grid minor;
end

function [Sig02_m E_1000_m Eps_m LamdaD_m SigRsurf_m] = fLoadDataFromFile(fileName)
    %% LAY SO LIEU TU FILE
    fileID = fopen(fileName, 'r');
    A=fscanf(fileID,'%f %f',[5 Inf]);
    fclose(fileID);
    B=A';

    Sig02_m=B(:,1);
    E_1000_m=B(:,2);
    Eps_m=B(:,3);
    LamdaD_m=B(:,4);
    SigRsurf_m=B(:,5);
end

%% 1. Eps (%)
function Eps = fEps(X)
    a1=-0.10285;
    a2=0.97991;

%     X=1000*Sig02/E;
    Eps=a1*X+a2;
end
%% 2. LamdaD (%)
function LamdaD = fLamdaD(X)
    a1=-0.10016;
    a2=0.000133539;

%     X=1000*Sig02/E;
    LamdaD=a1*X+a2;
end

%% 3. 1000*SigRsurf/E (%)
function SigRsurf_E = fSigRsurf_E(X)
    a1=0.21392;
    a2=-0.03701;

%     X=1000*Sig02/E;
    SigRsurf_E=a1*X+a2;
end