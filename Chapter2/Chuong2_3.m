clc; clear all; close all;

% %% 1. VE DO THI UNG SUAT KEO VA DUONG CONG CHAY
% Etb=21000; % MPa
% SigmaTtb=519; % MPa
% SigmaBtb=649; % MPa
% Del5tb=20.5; % %
% EpsP_Lim=0.02;
% 
% EpsT=SigmaTtb/Etb;
% Eps1=[0:0.1:1]*EpsT;
% Sigma1=Etb*Eps1;
% 
% EpsP=[0:0.0001:0.5];
% for i=1:length(EpsP)
%     SigmaS(i) = fSigmaS(EpsP(i),Etb,SigmaTtb,SigmaBtb,Del5tb,EpsP_Lim);
%     Eps2(i)=EpsP(i)+SigmaS(i)/Etb;
% end
% 
% Eps=[Eps1 Eps2];
% Sigma=[Sigma1 SigmaS];
% 
% xlabels{1} = 'Strain \epsilon [-]';
% xlabels{2} = 'Plastic strain \epsilon_{p} [-]';
% ylabels{1} = 'Stress \sigma [MPa]';
% ylabels{2} = 'Flow stress \sigma_{s} [MPa]';
% figure(1);
% [ax,hlT,hlS] = plotxx(Eps,Sigma,EpsP,SigmaS,xlabels,ylabels);

%% 2. KHAO SAT SU KHONG DONG DEU CUA PHAN BO US-BD
% function Main()
    Etb=210000; % MPa
    dE=5000;
    SigmaTtb=519; % MPa
    dSigmaT=59;
    SigmaBtb=649; % MPa
    dSigmaB=59;
    Del5tb=20; % %
    dDel5=1;
    Miutb=0.27;
    dMiu=0.1;

    n=12;
    Rtr=680;
    Rop=666;
    u=20;

    nMau=10^4; % So luong mau
    nBin=5; % So cot tren do thi hist

    kE=1;
    kSigmaT=2;
    kSigmaB=2;
    kDel5=0;
    kMiu=2;
    
    AlpE=0;
    AlpSigmaT=0;
    AlpSigmaB=0;
    AlpDel5=0;
    AlpMiu=0;

    kSigmaLim=1.07;
    kEpsLim=20;

    [E_1 SigmaT_1 SigmaB_1 Del5_1 Miu_1 kSigma_1 kEps_1...
     rE_1 rSigmaT_1 rSigmaB_1 rDel5_1 rMiu_1...
     D11_1 D12_1 D21_1 D22_1] = fMain(nMau,...
                                      Etb,dE,SigmaTtb,dSigmaT,SigmaBtb,dSigmaB,Del5tb,dDel5,Miutb,dMiu,...
                                      n,Rtr,Rop,u,...
                                      kE,AlpE,kSigmaT,AlpSigmaT,kSigmaB,AlpSigmaB,kDel5,AlpDel5,kMiu,AlpMiu,...
                                      kSigmaLim,kEpsLim);
    AlpE=0.4;
    AlpSigmaT=0.4;
    AlpSigmaB=0;
    AlpDel5=-0.2;
    AlpMiu=-0.4;
    
    [E_2 SigmaT_2 SigmaB_2 Del5_2 Miu_2 kSigma_2 kEps_2...
     rE_2 rSigmaT_2 rSigmaB_2 rDel5_2 rMiu_2...
     D11_2 D12_2 D21_2 D22_2] = fMain(nMau,...
                                      Etb,dE,SigmaTtb,dSigmaT,SigmaBtb,dSigmaB,Del5tb,dDel5,Miutb,dMiu,...
                                      n,Rtr,Rop,u,...
                                      kE,AlpE,kSigmaT,AlpSigmaT,kSigmaB,AlpSigmaB,kDel5,AlpDel5,kMiu,AlpMiu,...
                                      kSigmaLim,kEpsLim);

    %% VE DO THI
    figure(1);
    subplot(5,2,1);
    hist(E_1/10^3,nBin);
    % histogram(E,nBin,'BinWidth',(Emax-Emin)/nBin);
    xlabel('E [GPa]');
    ylabel('Count');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    subplot(5,2,2);
    hist(E_2/10^3,nBin);
    % histogram(E,nBin,'BinWidth',(Emax-Emin)/nBin);
    xlabel('E [GPa]');
    ylabel('Count');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;

    subplot(5,2,3);
    hist(SigmaT_1,nBin);
    % histogram(SigmaT,nBin,'BinWidth',(SigmaTmax-SigmaTmin)/nBin);
    xlabel('\sigma_{0.2} [MPa]');
    ylabel('Count');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    subplot(5,2,4);
    hist(SigmaT_2,nBin);
    % histogram(SigmaT,nBin,'BinWidth',(SigmaTmax-SigmaTmin)/nBin);
    xlabel('\sigma_{0.2} [MPa]');
    ylabel('Count');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;

    subplot(5,2,5);
    hist(SigmaB_1,nBin);
    % histogram(SigmaB,nBin,'BinWidth',(SigmaBmax-SigmaBmin)/nBin);
    xlabel('\sigma_{b} [MPa]');
    ylabel('Count');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    subplot(5,2,6);
    hist(SigmaB_2,nBin);
    % histogram(SigmaB,nBin,'BinWidth',(SigmaBmax-SigmaBmin)/nBin);
    xlabel('\sigma_{b} [MPa]');
    ylabel('Count');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;

    subplot(5,2,7);
    hist(Del5_1,nBin);
    % histogram(Del5,nBin,'BinWidth',(Del5max-Del5min)/nBin);
    xlabel('\delta_{5} [%]');
    ylabel('Count');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    subplot(5,2,8);
    hist(Del5_2,nBin);
    % histogram(Del5,nBin,'BinWidth',(Del5max-Del5min)/nBin);
    xlabel('\delta_{5} [%]');
    ylabel('Count');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;

    subplot(5,2,9);
    hist(Miu_1,nBin);
    % histogram(Miu,nBin,'BinWidth',(Miumax-Miumin)/nBin);
    xlabel('\mu [-]');
    ylabel('Count');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    subplot(5,2,10);
    hist(Miu_2,nBin);
    % histogram(Miu,nBin,'BinWidth',(Miumax-Miumin)/nBin);
    xlabel('\mu [-]');
    ylabel('Count');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;

    figure(2);
    subplot(2,2,1);
    hist(kSigma_1,nBin);
    xlabel('k_{\sigma} [-]');
    ylabel('Count');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    subplot(2,2,2);
    hist(kSigma_2,nBin);
    xlabel('k_{\sigma} [-]');
    ylabel('Count');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;

    subplot(2,2,3);
    hist(kEps_1,nBin);
    xlabel('k_{\epsilon} [-]');
    ylabel('Count');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    subplot(2,2,4);
    hist(kEps_2,nBin);
    xlabel('k_{\epsilon} [-]');
    ylabel('Count');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
% end

function [E SigmaT SigmaB Del5 Miu kSigma kEps...
          rE rSigmaT rSigmaB rDel5 rMiu...
          D11 D12 D21 D22] = fMain(nMau,...
                                   Etb,dE,SigmaTtb,dSigmaT,SigmaBtb,dSigmaB,Del5tb,dDel5,Miutb,dMiu,...
                                   n,Rtr,Rop,u,...
                                   kE,AlpE,kSigmaT,AlpSigmaT,kSigmaB,AlpSigmaB,kDel5,AlpDel5,kMiu,AlpMiu,...
                                   kSigmaLim,kEpsLim)

    Emin=Etb-dE;
    Emax=Etb+dE;
    SigmaTmin=SigmaTtb-dSigmaT;
    SigmaTmax=SigmaTtb+dSigmaT;
    SigmaBmin=SigmaBtb-dSigmaB;
    SigmaBmax=SigmaBtb+dSigmaB;
    Del5min=Del5tb-dDel5;
    Del5max=Del5tb+dDel5;
    Miumin=Miutb-dMiu;
    Miumax=Miutb+dMiu;
   
    Tetasao=fTetasao(n,Rtr,Rop,u);
 
    % Dij: i-sigma, j-eps; i,j=1-Tot; i,j=2-Khong tot
    D11=0;
    D12=0;
    D21=0;
    D22=0;

    i=0;
    while i<nMau
        SigmaTtmp=fPB(kSigmaT,SigmaTmin,SigmaTmax,AlpSigmaT);
        SigmaBtmp=fPB(kSigmaB,SigmaBmin,SigmaBmax,AlpSigmaB);
        kSigmaT_B=SigmaTtmp/SigmaBtmp;
        if kSigmaT_B<=0.9
            i=i+1;
            E(i)=fPB(kE,Emin,Emax,AlpE);
            SigmaT(i)=SigmaTtmp;
            SigmaB(i)=SigmaBtmp;
            Del5(i)=fPB(kDel5,Del5min,Del5max,AlpDel5);
            Miu(i)=fPB(kMiu,Miumin,Miumax,AlpMiu);

            kSigma(i)=exp(Tetasao*Miu(i));

            B=(log(SigmaB(i))-log(SigmaT(i)))/(log(Del5(i)/100-SigmaB(i)/E(i))-log(0.002));
            kEps(i)=(kSigma(i))^(1/B);

            if kSigma(i)<=kSigmaLim && kEps(i)<=kEpsLim
                D11=D11+1;
            elseif kSigma(i)<=kSigmaLim && kEps(i)>kEpsLim
                D12=D12+1;
            elseif kSigma(i)>kSigmaLim && kEps(i)<=kEpsLim
                D21=D21+1;
            else
                D22=D22+1;
            end
        end
    end
    
    D11=D11
    D12=D12
    D21=D21
    D22=D22
    
    kSigmaMax=max(kSigma);
    kSigmaMin=min(kSigma);
    kEpsMax=max(kEps);
    kEpsMin=min(kEps);

    %% TINH CAC HE SO TUONG QUAN
    rE = corr2(E,kEps)
    rSigmaT = corr2(SigmaT,kEps)
    rSigmaB = corr2(SigmaB,kEps)
    rDel5 = corr2(Del5,kEps)
    rMiu = corr2(Miu,kEps)
end


%% CAC HAM TOAN LIEN QUAN
function SigmaS = fSigmaS(EpsP,E,SigmaT,SigmaB,Del5,EpsP_Lim) % EpsP_Lim: Diem chay deo quy uoc (vi du: 0.002 - tuong ung 0.2%)
    EpsB=Del5/100;%log(1+Del5/100);
    EpsP_B=EpsB-SigmaB/E;
    N=log(SigmaB/SigmaT)/log(EpsP_B/EpsP_Lim);
    A=SigmaT/(EpsP_Lim)^N;
    
    if EpsP<=EpsP_Lim
        SigmaS=SigmaT;
    elseif EpsP<=EpsP_B
        SigmaS=A*EpsP^N;
    else
        SigmaS=SigmaB;
    end
end

function Tetasao = fTetasao(n,Rtr,Rop,u)
    Alp=2*pi/n;
    cosTetasao=(u+Rtr+Rop*(cos(Alp/2)-1))/sqrt((u+Rtr-Rop)^2+Rop^2+2*Rop*(u+Rtr-Rop)*cos(Alp/2));
    Tetasao=acos(cosTetasao);
end

%% HAM TAO PHAN BO
% Lua chon loai phan bo: 0-deu, 1-Simson, con lai-Chuan
% Alp: he so lech truc
function PB = fPB(kChon,xMin,xMax,Alp)
    k=Alp*(xMax-xMin)/2;
    PBtmp=xMax+1;
    if kChon==0
        while (PBtmp>xMax || PBtmp<xMin)
            a1=rand;
            PBtmp=(xMax-xMin)*a1+xMin+k;
        end
        PB=PBtmp;
    elseif kChon==1
        while (PBtmp>xMax || PBtmp<xMin)
            a=rand;
            b=rand;
            PBtmp=(xMax-xMin)*(a+b)/2+xMin+k;
        end
        PB=PBtmp;
    else
        while (PBtmp>xMax || PBtmp<xMin)
            s=0;
            for i=1:12
                s=s+rand;
            end
            s=s-6;
            PBtmp=(xMax+xMin)/2+(xMax-xMin)*s/6+k;
        end
        PB=PBtmp;
    end
end

function [ax,hl1,hl2] = plotxx(x1,y1,x2,y2,xlabels,ylabels)
    if nargin < 4
        error('Not enough input arguments')
    elseif nargin==4
        %Use empty strings for the xlabels
        xlabels{1}=' '; xlabels{2}=' '; ylabels{1}=' '; ylabels{2}=' ';
    elseif nargin==5
        %Use empty strings for the ylabel
        ylabels{1}=' '; ylabels{2}=' ';
    elseif nargin > 6
        error('Too many input arguments');
    end
    if length(ylabels) == 1
        ylabels{2} = ' ';
    end
    if ~iscellstr(xlabels) 
        error('Input xlabels must be a cell array');
    elseif ~iscellstr(ylabels) 
        error('Input ylabels must be a cell array');
    end
    hl1=line(x1,y1,'Color','k','LineWidth',2);
    ax(1)=gca;
    set(ax(1),'Position',[0.12 0.12 0.75 0.70]);
    ax(1).XAxis.FontSize = 12;
    ax(1).YAxis.FontSize = 12;
    set(ax(1),'XColor','k','YColor','k');
    ax(2)=axes('Position',get(ax(1),'Position'),...
               'XAxisLocation','top',...
               'YAxisLocation','right',...
               'Color','none',...
               'XColor','b','YColor','b');
    ax(2).XAxis.FontSize = 12;
    ax(2).YAxis.FontSize = 12;
    set(ax,'box','off');
    hl2=line(x2,y2,'Color','b','Parent',ax(2),'LineWidth',2);
    %label the two x-axes
    set(get(ax(1),'xlabel'),'string',xlabels{1});
    set(get(ax(2),'xlabel'),'string',xlabels{2});
    set(get(ax(1),'ylabel'),'string',ylabels{1});
    set(get(ax(2),'ylabel'),'string',ylabels{2});
    grid on;
    grid minor;
end