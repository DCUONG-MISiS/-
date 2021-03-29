function FuncStart
    clc; clear all; close all;
   
    tic;
    [Dn0MN OV0MN ...                                                                  % Hang so ban dau
         Dn1420MN OV1420MN Eps1420MN LamdaDMN SigT1420maxMN SigI1420maxMN ...% Hang so khi Dn=1420
         AlpD0MN D0MN AlpD1420MN D1420MN ...                                                   % Theo Teta 0-180
         Alp1420_1MN SigI1420MN Alp1420_2MN EpsT1420MN SigT1420MN  ...                                        % Theo Teta 0-360
         EpsMN DnMN OVMN ...                                 % Theo Eps
         ] =LoadData("MN5.txt");

    %% GRAPHICs
    % 1. Phan bo Duong kinh (D) khi Dn=1420mm
    figure(1);
    plot(AlpD0MN,D0MN,'b','LineWidth',1.5);
    hold on;
    plot(AlpD1420MN,D1420MN,'r','LineWidth',1.5);
    xticks([0:1/6:1]*180);
    xlim([0 180]);
    xlabel('\theta [град.]');
    ylabel('D [мм]');
    legend({'До экспандирования', 'После экспандирования'},'Box','off','FontSize',12);
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    % 2. Phan bo Cuong do Ung suat (SigI) khi Dn=1420mm - truoc khi do tai!
    figure(2);
    subplot(2,1,1);
    plot(Alp1420_1MN,SigI1420MN,'b','LineWidth',1.5);
    xticks([0:1/12:1]*360);
    xlim([0 360]);
    xlabel('\theta [град.]');
    ylabel('\sigma_{isurf} [МПа]');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    subplot(2,1,2);
    plot(Alp1420_2MN,SigT1420MN,'b','LineWidth',1.5);
    xticks([0:1/12:1]*360);
    xlim([0 360]);
    xlabel('\theta [град.]');
    ylabel('\sigma_{\thetarsurf} [МПа]');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    % 3. Phan bo Bien dang du (EpsT) khi Dn=1420mm - Sau khi do tai!
    figure(3);
    plot(Alp1420_2MN,EpsT1420MN,'b','LineWidth',1.5);
    xticks([0:1/12:1]*360);
    xlim([0 360]);
    xlabel('\theta [град.]');
    ylabel('\epsilon_{t} [%]');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;

    figure(4);
    yyaxis left
    plot(EpsMN,DnMN,'k','LineWidth',1.5);
    yyaxis right
    plot(EpsMN,OVMN,'b','LineWidth',1.5);
    xlabel('\epsilon_{ex} [%]');
    yyaxis left;
    ylabel('D_{n} [мм]');
    yyaxis right;
    ylabel('\Delta [мм]');
    ax = gca;
    ax.FontSize = 12;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'b';
    grid on;
    grid minor;

    [Dn0QT OV0QT ...                                                                  % Hang so ban dau
         Dn1420QT OV1420QT Eps1420QT LamdaDQT SigT1420maxQT SigI1420maxQT ...% Hang so khi Dn=1420
         AlpD0QT D0QT AlpD1420QT D1420QT ...                                                   % Theo Teta 0-180
         Alp1420_1QT SigI1420QT Alp1420_2QT EpsT1420QT SigT1420QT  ...                                        % Theo Teta 0-360
         EpsQT DnQT OVQT ...                                 % Theo Eps
         ] =LoadData("QT8.txt");

    %% GRAPHICs
    % 1. Phan bo Duong kinh (D) khi Dn=1420mm
    figure(5);
    plot(AlpD0QT,D0QT,'b','LineWidth',1.5);
    hold on;
    plot(AlpD1420QT,D1420QT,'r','LineWidth',1.5);
    xticks([0:1/6:1]*180);
    xlim([0 180]);
    xlabel('\theta [град.]');
    ylabel('D [мм]');
    legend({'До экспандирования', 'После экспандирования'},'Box','off','FontSize',12);
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    % 2. Phan bo Cuong do Ung suat (SigI) khi Dn=1420mm - truoc khi do tai!
    figure(6);
    subplot(2,1,1);
    plot(Alp1420_1QT,SigI1420QT,'b','LineWidth',1.5);
    xticks([0:1/12:1]*360);
    xlim([0 360]);
    xlabel('\theta [град.]');
    ylabel('\sigma_{isurf} [МПа]');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    subplot(2,1,2);
    plot(Alp1420_2QT,SigT1420QT,'b','LineWidth',1.5);
    xticks([0:1/12:1]*360);
    xlim([0 360]);
    xlabel('\theta [град.]');
    ylabel('\sigma_{\thetarsurf} [МПа]');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    % 3. Phan bo Bien dang du (EpsT) khi Dn=1420mm - Sau khi do tai!
    figure(7);
    plot(Alp1420_2QT,EpsT1420QT,'b','LineWidth',1.5);
    xticks([0:1/12:1]*360);
    xlim([0 360]);
    xlabel('\theta [град.]');
    ylabel('\epsilon_{t} [%]');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;

    figure(8);
    yyaxis left
    plot(EpsQT,DnQT,'k','LineWidth',1.5);
    yyaxis right
    plot(EpsQT,OVQT,'b','LineWidth',1.5);
    xlabel('\epsilon_{ex} [%]');
    yyaxis left;
    ylabel('D_{n} [мм]');
    yyaxis right;
    ylabel('\Delta [мм]');
    ax = gca;
    ax.FontSize = 12;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'b';
    grid on;
    grid minor;
    
    figure(9);
    subplot(2,1,1);
    plot(Alp1420_1MN,SigI1420MN,'b','LineWidth',1.5);
    xticks([0:1/12:1]*360);
    xlim([0 360]);
    xlabel('\theta [град.]');
    ylabel('\sigma_{isurf} [МПа]');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    subplot(2,1,2);
    plot(Alp1420_1QT,SigI1420QT,'b','LineWidth',1.5);
    xticks([0:1/12:1]*360);
    xlim([0 360]);
    xlabel('\theta [град.]');
    ylabel('\sigma_{isurf} [МПа]');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    figure(10);
    subplot(2,1,1);
    plot(Alp1420_2MN,SigT1420MN,'b','LineWidth',1.5);
    xticks([0:1/12:1]*360);
    xlim([0 360]);
    xlabel('\theta [град.]');
    ylabel('\sigma_{\thetarsurf} [МПа]');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    subplot(2,1,2);
    plot(Alp1420_2QT,SigT1420QT,'b','LineWidth',1.5);
    xticks([0:1/12:1]*360);
    xlim([0 360]);
    xlabel('\theta [град.]');
    ylabel('\sigma_{\thetarsurf} [МПа]');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    grid on;
    grid minor;
    
    toc;
end

%% Ham lay so lieu tu mo phong
function [Dn0 OV0 ...                                                                  % Hang so ban dau
         Dn1420 OV1420 Eps1420 LamdaD SigT1420max SigI1420max ...% Hang so khi Dn=1420
         AlpD0 D0 AlpD1420 D1420 ...                                                   % Theo Teta 0-180
         Alp1420_1 SigI1420 Alp1420_2 EpsT1420 SigT1420  ...                                        % Theo Teta 0-360
         Eps Dn OV ...                                 % Theo Eps
         ] =LoadData(FileName)

    %% LAY SO LIEU TU FILE
    fileID = fopen(FileName, 'r');
    A=fscanf(fileID,'%f %f',[11 Inf]);
    fclose(fileID);
    B=A';
    
    RecTmp=B(:,1);
    XTmp=B(:,2);
    ZTmp=B(:,3);
    AlpTmp=B(:,4);
    RTmp=B(:,5);
    EpsPTmp=B(:,6);
    EpsRTmp=B(:,7);
    EpsTTmp=B(:,8);
    SigITmp=B(:,9);
    SigRTmp=B(:,10);
    SigTTmp=B(:,11);
    
    cRec=hist(RecTmp,1:max(RecTmp));
    ncRec=length(cRec);
    k=0;
    for i=1:ncRec
        for j=1:cRec(i)
            X{i}(j)=XTmp(j+k);
            Z{i}(j)=ZTmp(j+k);
            Alp{i}(j)=AlpTmp(j+k);
            R{i}(j)=RTmp(j+k);
            EpsP{i}(j)=EpsPTmp(j+k);
            EpsR{i}(j)=EpsRTmp(j+k);
            EpsT{i}(j)=EpsTTmp(j+k);
            SigI{i}(j)=SigITmp(j+k);
            SigR{i}(j)=SigRTmp(j+k);
            SigT{i}(j)=SigTTmp(j+k);
        end
        k=k+cRec(i);
    end

    %% SAP XEP THEO THU TU TANG DAN CUA GOC ALP
    for i=1:ncRec
        [tmpAlp{i} indAlp{i}]=sort(Alp{i});
        Nr=cRec(i);
        for j=1:Nr
            tmpX{i}(j)=X{i}(indAlp{i}(j));
            tmpZ{i}(j)=Z{i}(indAlp{i}(j));
            tmpR{i}(j)=R{i}(indAlp{i}(j));
            tmpEpsP{i}(j)=EpsP{i}(indAlp{i}(j));
            tmpEpsR{i}(j)=EpsR{i}(indAlp{i}(j));
            tmpEpsT{i}(j)=EpsT{i}(indAlp{i}(j));
            tmpSigI{i}(j)=SigI{i}(indAlp{i}(j));
            tmpSigR{i}(j)=SigR{i}(indAlp{i}(j));
            tmpSigT{i}(j)=SigT{i}(indAlp{i}(j));
        end
        Alp{i}=tmpAlp{i};
        X{i}=tmpX{i};
        Z{i}=tmpZ{i};
        R{i}=tmpR{i};
        EpsP{i}=100*tmpEpsP{i};
        EpsR{i}=100*tmpEpsR{i};
        EpsT{i}=100*tmpEpsT{i};
        SigI{i}=tmpSigI{i};
        SigR{i}=tmpSigR{i};
        SigT{i}=tmpSigT{i};
    end
    
    for i=1:ncRec
        AlpDSample=[0:0.1:180];
        nAlpD=length(AlpDSample);
        for j=1:nAlpD
            AlpD{i}(j)=AlpDSample(j);
            D{i}(j) = fD(Alp{i},R{i},AlpDSample(j));
        end
    end
    
    %% TINH CHUYEN VI U, BAN KINH Rmax, Rmin, CHU VI CV va DO TRON DTr
    for i=1:ncRec
        Dmax=max(D{i});
        Dmin=min(D{i});
        OV(i)=Dmax-Dmin;
        % Tinh chu vi
        Nr=length(Alp{i});
        CV(i)=0;
        for j=1:Nr-1
            CV(i)=CV(i)+sqrt((X{i}(j+1)-X{i}(j))^2+(Z{i}(j+1)-Z{i}(j))^2);
        end
        CV(i)=CV(i)+sqrt((X{i}(Nr)-X{i}(1))^2+(Z{i}(Nr)-Z{i}(1))^2);
        Eps(i)=100*(CV(i)-CV(1))/CV(1);
        
        Dn(i)=CV(i)/pi;
    end
    
    %% TINH CAC HANG SO
    Dn0=Dn(1);
    OV0=OV(1);

    Dn1420=Dn(ncRec);
    OV1420=OV(ncRec);
    Eps1420=Eps(ncRec);
    
    LamdaD=(Dn(ncRec)-Dn(ncRec-1))*100/Dn(ncRec-1);

%     EpsT1420min=min(EpsT{ncRec});
%     EpsT1420max=max(EpsT{ncRec});
%     
%     SigI1420min=min(SigI{ncRec-1});
    SigI1420max=max(SigI{ncRec-1});
%     SigT1420min=min(SigT{ncRec});
    SigT1420max=max(SigT{ncRec});
    
%     mSigI1420 = fMsig(Alp{ncRec-1}, SigI{ncRec-1});
    
    AlpD0=AlpD{1};
    D0=D{1};
    AlpD1420=AlpD{ncRec};
    D1420=D{ncRec};
    
    Alp1420_1=Alp{ncRec-1};
    SigI1420=SigI{ncRec-1};
    
    Alp1420_2=Alp{ncRec};
    EpsT1420=EpsT{ncRec};
    SigT1420=SigT{ncRec};
end

function [SigItb mSigI Eta] =LoadDataLastRecord(FileName)
    
    global Sig02

    %% LAY SO LIEU TU FILE
    fileID = fopen(FileName, 'r');
    A=fscanf(fileID,'%f',[1 Inf]);
    fclose(fileID);
    B=A';
    
    SigI=B(:,1);
    SigItb=mean(SigI);
    
    dem=0;
    N=length(SigI);
    for i=1:N
        if SigI(i)<=Sig02
            dem=dem+1;
        end
    end
    SigImax=max(SigI);
    mSigI=SigImax/SigItb;
    Eta=dem/N;
end

function rSigR =LoadData_rSigR(FileName)
    %% LAY SO LIEU TU FILE
    fileID = fopen(FileName, 'r');
    A=fscanf(fileID,'%f',[1 Inf]);
    fclose(fileID);
    B=A';
    
    SigR=B(:,1);
    
    dem=0;
    N=length(SigR);
    for i=1:N
        if SigR(i)>0
            dem=dem+1;
        end
    end
    rSigR=dem*100/N;
end

function mSig_Eps=fmSig_Eps(X)
    n=length(X);
    dem=0;
    for i=1:n
        if X(i)>0
            dem=dem+1;
        end
    end
    mSig_Eps=dem*100/n;
end

function vRD = fRD(Alp,R,AlpD)
    n=length(Alp);
    if AlpD<Alp(1)
        vRD=R(1)+(R(2)-R(1))*(AlpD-Alp(1))/(Alp(2)-Alp(1));
    elseif AlpD>Alp(n)
        vRD=R(n-1)+(R(n)-R(n-1))*(AlpD-Alp(n-1))/(Alp(n)-Alp(n-1));
    else
        vRD=interp1(Alp,R,AlpD);
    end
end

function vD = fD(Alp,R,AlpD)
    k=pi/180;
    RD1 = fRD(Alp,R,AlpD);
    X1=RD1*cos(k*AlpD);
    Z1=RD1*sin(k*AlpD);
    RD2 = fRD(Alp,R,AlpD+180);
    X2=RD2*cos(k*(AlpD+180));
    Z2=RD2*sin(k*(AlpD+180));
    vD=sqrt((X1-X2)^2+(Z1-Z2)^2);
end

function mSigI = fMsig(Alp, SigI)
    n=length(Alp);
    for i=1:n+2
        if i==1
            Alptmp(i)=0;
            SigItmp(i)=SigI(1);
        elseif i==n+2
            Alptmp(i)=360;
            SigItmp(i)=SigI(n);
        else
            Alptmp(i)=Alp(i-1);
            SigItmp(i)=SigI(i-1);
        end
    end
    s=0;
    for i=1:n+1
        s=s+(SigItmp(i+1)+SigItmp(i))*(Alptmp(i+1)-Alptmp(i))/2;
    end
    SigImax=max(SigI);
    SigItb=s/360;
    if SigItb==0
        mSigI=0;
    else
        mSigI=SigImax/SigItb;
    end
end

function plot3y(x1,y1,x2,y2,x3,y3,xlbel,ylbels)

    if nargin==7
       %Use empty strings for the ylabels
       ylabels{1}=' '; ylabels{2}=' '; ylabels{3}=' ';
    elseif nargin > 8
       error('Too many input arguments')
    elseif nargin < 7
       error('Not enough input arguments')
    end

    yyaxis left
    plot(x1,y1,'k','LineWidth',1.5);
    yyaxis right
    plot(x2,y2,'b','LineWidth',1.5);
    xlabel(xlbel);
    yyaxis left;
    ylabel(ylbels{1});
    yyaxis right;
    ylabel(ylbels{2});
    ax = gca;
    ax.FontSize = 12;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'b';
    grid on;
    grid minor;
    
    pos =get(ax,'position');
    kW=5;
    dW=pos(3)/kW;
    %Reduce width of the two axes generated by plotyy 
    pos(3) = pos(3)-dW;
    set(ax,'position',pos);  
    %Determine the position of the third axes
    pos3=pos+[0 0 dW 0];
    %Determine the proper x-limits for the third axes
    limx1=get(ax,'xlim');
    limx3=[limx1(1) limx1(1)+kW*(limx1(2)-limx1(1))/(kW-1)];

    ax(3)=axes('Position',pos3,'box','off',...
       'Color','none','XColor','k','YColor','r',...   
       'xtick',[],'xlim',limx3,'yaxislocation','right');
    line(x3,y3,'Color','r','LineWidth',1.5,'Parent',ax(3));
    ax(3).XAxis.Visible = 'off';
    set(get(ax(3),'ylabel'),'string',ylbels{3})
    ax(3).YAxis.FontSize = 12;
end

function plot4y(x1,y1,x2,y2,x3,y3,x4,y4,xlbel,ylbels)

    if nargin==9
       %Use empty strings for the ylabels
       ylabels{1}=' '; ylabels{2}=' '; ylabels{3}=' '; ylabels{4}=' ';
    elseif nargin > 10
       error('Too many input arguments')
    elseif nargin < 9
       error('Not enough input arguments')
    end

    yyaxis left
    plot(x1,y1,'k','LineWidth',1.5);
    yyaxis right
    plot(x2,y2,'b','LineWidth',1.5);
    xlabel(xlbel);
    yyaxis left;
    ylabel(ylbels{1});
    yyaxis right;
    ylabel(ylbels{2});
    ax = gca;
    ax.FontSize = 12;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'b';
    grid on;
    grid minor;
    
    pos =get(ax,'position');
    kW=5;
    dW=pos(3)/kW;
    %Reduce width of the two axes generated by plotyy 
    pos(3) = pos(3)-2*dW;
    set(ax,'position',pos);  
    %Determine the position of the third axes
    pos3=pos+[0 0 dW 0];
    pos4=pos+[0 0 2*dW 0];
    %Determine the proper x-limits for the third axes
    limx1=get(ax,'xlim');
    limx3=[limx1(1) limx1(1)+(kW-1)*(limx1(2)-limx1(1))/(kW-2)];
    limx4=[limx1(1) limx1(1)+kW*(limx1(2)-limx1(1))/(kW-2)];

    ax(3)=axes('Position',pos3,'box','off',...
       'Color','none','XColor','k','YColor','r',...   
       'xtick',[],'xlim',limx3,'yaxislocation','right');
    line(x3,y3,'Color','r','LineWidth',1.5,'Parent',ax(3));
    ax(3).XAxis.Visible = 'off';
    set(get(ax(3),'ylabel'),'string',ylbels{3})
    ax(3).YAxis.FontSize = 12;
    
    ax(4)=axes('Position',pos4,'box','off',...
       'Color','none','XColor','k','YColor','m',...   
       'xtick',[],'xlim',limx4,'yaxislocation','right');
    line(x4,y4,'Color','m','LineWidth',1.5,'Parent',ax(4));
    ax(4).XAxis.Visible = 'off';
    set(get(ax(4),'ylabel'),'string',ylbels{4})
    ax(4).YAxis.FontSize = 12;
end