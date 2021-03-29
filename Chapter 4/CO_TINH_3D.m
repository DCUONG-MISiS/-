function Batdau()
    close all; clc; clear all;

    [Sig02_m E_m Eps_m LamdaD_m SigRsurf_m] = fLoadDataFromFile('Cotinh.txt');

    N=50;
    Sig02min=min(Sig02_m);
    Sig02max=max(Sig02_m);
    Sig02=[Sig02min:(Sig02max-Sig02min)/N:Sig02max];
    Emin=min(E_m);
    Emax=max(E_m);
    E=[Emin:(Emax-Emin)/N:Emax];
    
    for i=1:N+1
        for j=1:N+1
            Eps(i,j) = fEps(Sig02(i),E(j));
            LamdaD(i,j) = fLamdaD(Sig02(i),E(j));
            SigRsurf(i,j) = fSigRsurf(Sig02(i),E(j));
        end
    end

    %% 1. ///
    figure(1);
    scatter3(E_m,Sig02_m,Eps_m,60,'^','filled');
    hold on;
    surf(E,Sig02,Eps,'EdgeColor','none');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.XAxis.Color = 'b';
    ax.YAxis.FontSize = 12;
    ax.YAxis.Color = 'k';
    ax.ZAxis.FontSize = 12;
    ax.ZAxis.Color = 'r';
    xlabel('E [ГПа]','FontSize',16);
    ylabel('\sigma_{0,2} [МПа]','FontSize',16);
    zlabel('\epsilon [%]','FontSize',16);
    legend({'Данные моделирования','Аппроксимация'},'FontSize',12);
    grid on;
    grid minor; 
    
    %% 2./////
    figure(2);
    scatter3(E_m,Sig02_m,LamdaD_m,60,'^','filled');
    hold on;
    surf(E,Sig02,LamdaD,'EdgeColor','none');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.XAxis.Color = 'b';
    ax.YAxis.FontSize = 12;
    ax.YAxis.Color = 'k';
    ax.ZAxis.FontSize = 12;
    ax.ZAxis.Color = 'r';
    xlabel('E [ГПа]','FontSize',16);
    ylabel('\sigma_{0,2} [МПа]','FontSize',16);
    zlabel('\lambda_{D} [%]','FontSize',16);
    legend({'Данные моделирования','Аппроксимация'},'FontSize',12);
    grid on;
    grid minor;

    %% 3./////
    figure(3);
    scatter3(E_m,Sig02_m,SigRsurf_m,60,'^','filled');
    hold on;
    surf(E,Sig02,SigRsurf,'EdgeColor','none');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.XAxis.Color = 'b';
    ax.YAxis.FontSize = 12;
    ax.YAxis.Color = 'k';
    ax.ZAxis.FontSize = 12;
    ax.ZAxis.Color = 'r';
    xlabel('E [ГПа]','FontSize',16);
    ylabel('\sigma_{0,2} [МПа]','FontSize',16);
    zlabel('\sigma_{rsurf} [МПа]','FontSize',16);
    legend({'Данные моделирования','Аппроксимация'},'FontSize',12);
    grid on;
    grid minor;
end

function [Sig02_m E_m Eps_m LamdaD_m SigRsurf_m] = fLoadDataFromFile(fileName)
    %% LAY SO LIEU TU FILE
    fileID = fopen(fileName, 'r');
    A=fscanf(fileID,'%f %f',[5 Inf]);
    fclose(fileID);
    B=A';

    Sig02_m=B(:,1);
    E_m=B(:,2);
    Eps_m=B(:,3);
    LamdaD_m=B(:,4);
    SigRsurf_m=B(:,5);
end

%% 1. Eps (%)
function Eps = fEps(Sig02, E)
    a1=-0.10285;
    a2=0.97991;
    Eps=a1*Sig02/E+a2;
end
%% 2. LamdaD (%)
function LamdaD = fLamdaD(Sig02, E)
    a1=-0.10016;
    a2=0.000133539;

    LamdaD=a1*Sig02/E+a2;
end

%% 3. SigRsurf (MPa)
function SigRsurf = fSigRsurf(Sig02, E)
    a1=0.21392;
    a2=-0.03701;

    SigRsurf=a1*Sig02+a2*E;
end