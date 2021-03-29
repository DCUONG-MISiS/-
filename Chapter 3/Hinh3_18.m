% x1=[2.5:0.05:4.5];
% y3=interp1(x,y1,x1,'PCHIP');
% y4=interp1(x,y2,x1,'PCHIP');
% for i=1:length(x)
%     text(1.02*x(i),y1(i),num2str(y1(i)),'Color','k','FontSize',11);
%     text(1.02*x(i),y2(i),num2str(y2(i)),'Color','b','FontSize',11);
% end

clc; clear all; close all;

U=[0.5 1 1.5 2 2.5];

Emp1 = [1.86 2.81 4.51 6.74 7.80];
Del1=[7.91 5.80 -6.34 -4.26 -3.15];
Etn1 = Emp1./(1+Del1/100);
Cn1=188.4*(1+Etn1/100);

Emp2 = [1.59 2.65 4.35 5.94 7.59];
Del2=[8.99 7.43 3.19 -5.45 -2.83];
Etn2 = Emp2./(1+Del2/100);
Cn2=188.4*(1+Etn2/100);

Emp3 = [1.38 2.55 4.14 5.79 7.27];
Del3=[8.52 8.82 5.41 4.37 6.72];
Etn3 = Emp3./(1+Del3/100);
Cn3=188.4*(1+Etn3/100);

Emp4 = [1.59 2.92 4.35 5.84 7.64];
Del4=[7.32 8.03 7.84 -1.30 -5.56];
Etn4 = Emp4./(1+Del4/100);
Cn4=188.4*(1+Etn4/100);

Emp5 = [1.33 2.60 4.09 5.57 7.32];
Del5=[7.50 -6.04 -5.01 0.19 8.73];
Etn5 = Emp5./(1+Del5/100);
Cn5=188.4*(1+Etn5/100);

figure(1);
subplot(1,3,1);
plot(U,Emp1,'-kd','LineWidth',1.5);
hold on;
plot(U,Etn1,'-bs','LineWidth',1.5);
xticks([0.5:0.5:2.5]);
xlabel('Рабочий ход сегментов, U [мм]','FontWeight','bold');
ylabel('Суммарная деформация, \epsilon^{\Sigma}_{\theta} [%]','FontWeight','bold');
legend ({'Моделирование','Эксперимент'},'FontSize',12);
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
grid on;
grid minor;

subplot(1,3,2);
plot(U,Emp2,'-kd','LineWidth',1.5);
hold on;
plot(U,Etn2,'-bs','LineWidth',1.5);
xticks([0.5:0.5:2.5]);
xlabel('Рабочий ход сегментов, U [мм]','FontWeight','bold');
ylabel('Суммарная деформация, \epsilon^{\Sigma}_{\theta} [%]','FontWeight','bold');
legend ({'Моделирование','Эксперимент'},'FontSize',12);
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
grid on;
grid minor;

subplot(1,3,3);
plot(U,Emp3,'-kd','LineWidth',1.5);
hold on;
plot(U,Etn3,'-bs','LineWidth',1.5);
xticks([0.5:0.5:2.5]);
xlabel('Рабочий ход сегментов, U [мм]','FontWeight','bold');
ylabel('Суммарная деформация, \epsilon^{\Sigma}_{\theta} [%]','FontWeight','bold');
legend ({'Моделирование','Эксперимент'},'FontSize',12);
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
grid on;
grid minor;

figure(2);
subplot(1,2,1);
plot(U,Emp4,'-kd','LineWidth',1.5);
hold on;
plot(U,Etn4,'-bs','LineWidth',1.5);
xticks([0.5:0.5:2.5]);
xlabel('Рабочий ход сегментов, U [мм]','FontWeight','bold');
ylabel('Суммарная деформация, \epsilon^{\Sigma}_{\theta} [%]','FontWeight','bold');
legend ({'Моделирование','Эксперимент'},'FontSize',12);
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
grid on;
grid minor;

subplot(1,2,2);
plot(U,Emp5,'-kd','LineWidth',1.5);
hold on;
plot(U,Etn5,'-bs','LineWidth',1.5);
xticks([0.5:0.5:2.5]);
xlabel('Рабочий ход сегментов, U [мм]','FontWeight','bold');
ylabel('Суммарная деформация, \epsilon^{\Sigma}_{\theta} [%]','FontWeight','bold');
legend ({'Моделирование','Эксперимент'},'FontSize',12);
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
grid on;
grid minor;



