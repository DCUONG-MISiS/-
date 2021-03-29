% Long Name	Short Name	RGB Triplet
% blue	b	[0,0,1]
% black	k	[0,0,0]
% red	r	[1,0,0]
% green	g	[0,1,0]
% yellow	y	[1,1,0]
% cyan	c	[0,1,1]
% magenta	m	[1,0,1]
% white	w	[1,1,1]
% plot(X0, Y0, 'k', 'LineWidth', 3);
%  P=plot(Alp, Et, Beta,Ett);
%     P(1).LineWidth=1.5;
%     P(2).LineWidth=2;
% xticks([0:22.5:360]);
%     xlim([0 360]);
%     legend('11111111111111','11111111111112');
% legend boxoff;

function FuncStart
    clc; clear all; close all;
    
    global n Rc Rop u Miu K x s0 Rt Alp_2 Tetak sc Mc Sigc Epsc Tauc Pc nTeta1 nTeta2
    
    Miu=0.17;
    E=210*10^3;
    Sig02=460;
    Sigb=590;
    Del5=20;
    Epspb=Del5/100-Sigb/E;
    x=log(Sigb/Sig02)/log(Epspb/0.002);
    K=Sig02/0.002^x;
    
    n=12;
    u=10;
    Rc=670;
    Rop=670;
    s0=22;
    Rt=678;
    sc=s0-0.2;
    Mc=0;
    
    Alp_2=pi/n;
    cosTetak=(u+Rt+Rop*(cos(Alp_2)-1))/sqrt((u+Rt-Rop)^2+Rop^2+2*Rop*(u+Rt-Rop)*cos(Alp_2));
    Tetak=acos(cosTetak);
    
    Epsc=abs((sc-s0)/s0);
    Sigc=K*(2*Epsc/sqrt(3))^x;
    Tauc=Sigc*exp((Miu+1/Miu)*Tetak)*tan(Tetak);
    Pc=sc*(Sigc+Tauc/Miu);
    
    nTeta1=10;
    nTeta2=5;
    
    [Teta M] = fSolvM();
    [Teta1 P1 Teta2 P2] = fP_cal(Teta);
    nTeta=length(Teta);
    for i=1:nTeta
        s(i)=fS(Teta(i));
        Sig(i)=fSig(Teta(i));
        Tau(i)=fTau(Teta(i));
        Eps(i)=fEps(Teta(i));
        rSig(i) = frSig(Teta(i));
        rEps(i) = frEps(Teta(i));
    end
    
    fileID = fopen('SLLLL.txt', 'w');
    for i=1:nTeta
        if i<=nTeta1+1
            fprintf(fileID,'%d %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f\n',[i Teta(i) Sig(i) Tau(i) P1(i) s(i) Eps(i) M(i)/10^5]);
        else
            fprintf(fileID,'%d %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f\n',[i Teta(i) Sig(i) Tau(i) P2(i-nTeta1) s(i) Eps(i) M(i)/10^5]);
        end
    end
    fclose(fileID);
    
    TetaMax=max(Teta);
    
    figure(1);
    plot(Teta,Sig,'color','b','LineWidth',2);
    xticks([0:1/5:1]*TetaMax);
    xlim([0 TetaMax]);
    xlabel("Угол по 1/" + num2str(2*n) + " периметра, \theta [град.]",'FontWeight','bold');
    ylabel('Тангенциальное напряжение, \sigma [МПа]','FontWeight','bold');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    figure(2);
    plot(Teta,Tau,'color','b','LineWidth',2);
    xticks([0:1/5:1]*TetaMax);
    xlim([0 TetaMax]);
    xlabel("Угол по 1/" + num2str(2*n) + " периметра, \theta [град.]",'FontWeight','bold');
    ylabel('Касательное напряжение, \tau [MPa]','FontWeight','bold');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    figure(3);
    plot(Teta,s,'color','b','LineWidth',2);
    xticks([0:1/5:1]*TetaMax);
    xlim([0 TetaMax]);
    xlabel("Угол по 1/" + num2str(2*n) + " периметра, \theta [град.]",'FontWeight','bold');
    ylabel('Толщина стенки, s [мм]','FontWeight','bold');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    figure(4);
    plot(Teta,Eps,'color','b','LineWidth',2);
    xticks([0:1/5:1]*TetaMax);
    xlim([0 TetaMax]);
    xlabel("Угол по 1/" + num2str(2*n) + " периметра, \theta [град.]",'FontWeight','bold');
    ylabel('Тангенциальная деформация, \epsilon [-]','FontWeight','bold');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    figure(5);
    plot(Teta1,P1,Teta2,P2,'color','b','LineWidth',2);
    % Ve duong net dut
    hold on;
    xD=[Teta1(nTeta1+1) Teta2(1)];
    yD=[P1(nTeta1+1) P2(1)];
    plot(xD,yD,'b--','LineWidth',2);
    xticks([0:1/5:1]*TetaMax);
    xlim([0 TetaMax]);
    xlabel("Угол по 1/" + num2str(2*n) + " периметра, \theta [град.]",'FontWeight','bold');
    ylabel('Давление сегмента, p [МПа]','FontWeight','bold');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    figure(6);
    plot(Teta,M,'color','b','LineWidth',2);
    xticks([0:1/5:1]*TetaMax);
    xlim([0 TetaMax]);
    xlabel("Угол по 1/" + num2str(2*n) + " периметра, \theta [град.]",'FontWeight','bold');
    ylabel('Изгибающий момент, M [Н.м/м]','FontWeight','bold');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    figure(7);
    plot(Teta,rSig,'color','b','LineWidth',2);
    xticks([0:1/5:1]*TetaMax);
    xlim([0 TetaMax]);
    xlabel("Угол по 1/" + num2str(2*n) + " периметра, \theta [град.]",'FontWeight','bold');
    ylabel('Безразмерное напряжение, r_{\sigma} [-]','FontWeight','bold');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    figure(8);
    plot(Teta,rEps,'color','b','LineWidth',2);
    xticks([0:1/5:1]*TetaMax);
    xlim([0 TetaMax]);
    xlabel("Угол по 1/" + num2str(2*n) + " периметра, \theta [град.]",'FontWeight','bold');
    ylabel('Безразмерная деформация, r_{\epsilon} [-]','FontWeight','bold');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    n_=[10:2:20];
    Miu1=0.1;
    Miu2=0.2;
    Miu3=0.3;
    Miu4=0.4;
    for i=1:length(n_)
        kSig1(i)=fkSig(Miu1,n_(i));
        kSig2(i)=fkSig(Miu2,n_(i));
        kSig3(i)=fkSig(Miu3,n_(i));
        kSig4(i)=fkSig(Miu4,n_(i));
        
        kEps1(i)=fkEps(Miu1,n_(i));
        kEps2(i)=fkEps(Miu2,n_(i));
        kEps3(i)=fkEps(Miu3,n_(i));
        kEps4(i)=fkEps(Miu4,n_(i));
    end
    n_max=max(n_);
    n_min=min(n_);
    
    figure(9);
    plot(n_,kSig1,'color','b','LineWidth',2);
    hold on;
    plot(n_,kSig2,'color','k','LineWidth',2);
    hold on;
    plot(n_,kSig3,'color','y','LineWidth',2);
    hold on;
    plot(n_,kSig4,'color','r','LineWidth',2);
    xticks(n_);
    xlim([n_min n_max]);
    xlabel('Количество сегментов, n [шт]','FontWeight','bold');
    ylabel({'Степень неравномерности', 'напряжения, k_{\sigma} [-]'},'FontWeight','bold');
    legend({"\mu = "+num2str(Miu1),"\mu = "+num2str(Miu2),"\mu = "+num2str(Miu3),"\mu = "+num2str(Miu4)},'FontSize',12,'FontWeight','bold');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    figure(10);
    plot(n_,kEps1,'color','b','LineWidth',2);
    hold on;
    plot(n_,kEps2,'color','k','LineWidth',2);
    hold on;
    plot(n_,kEps3,'color','y','LineWidth',2);
    hold on;
    plot(n_,kEps4,'color','r','LineWidth',2);
    xticks(n_);
    xlim([n_min n_max]);
    xlabel('Количество сегментов, n [шт]','FontWeight','bold');
    ylabel({'Степень неравномерности', 'деформации, k_{\epsilon} [-]'},'FontWeight','bold');
    legend({"\mu = "+num2str(Miu1),"\mu = "+num2str(Miu2),"\mu = "+num2str(Miu3),"\mu = "+num2str(Miu4)},'FontSize',12,'FontWeight','bold');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
end

function fSig = fSig(Teta) % Teta: Do
    global Miu Sigc Tetak
    Tetatmp=pi*Teta/180;
    if Tetatmp<=Tetak
        fSig=Sigc*exp(Miu*Tetatmp);
    else
       fSig=Sigc*exp(Miu*Tetak)*cos(Tetatmp)/cos(Tetak);
    end
end

function fEps = fEps(Teta) % Teta: Do
    global K x
    Sig=fSig(Teta);
    fEps=(sqrt(3)/2)*(Sig/K)^(1/x);
end

function fP = fP(Teta) % Teta: Do
    global Miu Rc Tetak
    Tetatmp=pi*Teta/180;
    if Tetatmp<=Tetak
        s=fS(Teta);
        Sig=fSig(Teta);
        Tau=fTau(Teta);
        
        fP=s*(Sig+Tau/Miu)/Rc;
    else
       fP=0;
    end
end

function [Teta1 P1 Teta2 P2] = fP_cal(Teta) % Teta: Do
    global Tetak Alp_2 nTeta1 nTeta2
    dTeta1=Tetak/nTeta1;
    dTeta2=(Alp_2-Tetak)/nTeta2;
    for i=1:nTeta1+1
        Teta1(i)=(i-1)*dTeta1;
        P1(i)=fP(180*Teta1(i)/pi);
    end
    
    for i=1:nTeta2+1
        Teta2(i)=Tetak+(i-1)*dTeta2;
        P2(i)=0;
    end
    Teta1=180*Teta1/pi;
    Teta2=180*Teta2/pi;
end

function fS = fS(Teta) % Teta: Do
    global s0
    Eps=fEps(Teta);
    fS=s0*(1-Eps);
end

function fTau = fTau(Teta) % Teta: Do
    global Miu Tauc Tetak
    Tetatmp=pi*Teta/180;
    if Tetatmp<=Tetak
        fTau=Tauc*exp(-Tetatmp/Miu);
    else
       fTau=Tauc*exp(-Tetak/Miu)*sin(Tetatmp)/sin(Tetak);
    end
end

function dM = dM(Teta)
    global Miu Rc Tetak
    s=fS(Teta);
    P=fP(Teta);
    Tau=fTau(Teta);
    
    Tetatmp=pi*Teta/180;
    if Tetatmp<=Tetak
        dM=s*Rc*(Miu*P-2*Tau)/2;
    else
       dM=-s*Rc*Tau;
    end
end

function [Teta M] = fSolvM()
    global Mc Tetak Alp_2 nTeta1 nTeta2
    dTeta1=Tetak/nTeta1;
    dTeta2=(Alp_2-Tetak)/nTeta2;
    nTeta=nTeta1+nTeta2+1;
    for i=1:nTeta
        if i==1
            Teta(i)=0;
            M(i)=Mc;
        elseif i<=nTeta1+1
            Teta(i)=Teta(i-1)+dTeta1;
            M(i)=M(i-1)+dTeta1*dM(180*Teta(i-1)/pi);
        else
            Teta(i)=Teta(i-1)+dTeta2;
            M(i)=M(i-1)+dTeta2*dM(180*Teta(i-1)/pi);
        end
    end
    Teta=180*Teta/pi;
end


function rSig = frSig(Teta) % Teta: Do
    global Miu Tetak
    Tetatmp=pi*Teta/180;
    if Tetatmp<=Tetak
        rSig=exp(Miu*Tetatmp);
    else
       rSig=exp(Miu*Tetak)*cos(Tetatmp)/cos(Tetak);
    end
end

function kSig = fkSig(Miu_,n_)
    global Rt Rop u
    Alp_2=pi/n_;
    cosTetak=(u+Rt+Rop*(cos(Alp_2)-1))/sqrt((u+Rt-Rop)^2+Rop^2+2*Rop*(u+Rt-Rop)*cos(Alp_2));
    Tetak=acos(cosTetak);
    kSig=exp(Miu_*Tetak);
end

function rEps = frEps(Teta) % Teta: Do
    global x
    rSigtmp=frSig(Teta);
    rEps=(rSigtmp)^(1/x);
end

function kEps = fkEps(Miu_,n_)
    global x
    kSigtmp=fkSig(Miu_,n_);
    kEps=kSigtmp^(1/x);
end
