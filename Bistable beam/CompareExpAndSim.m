clc
clear all
close all

load('SingleBeamCompareDara.mat')

%%
plot(d_plot*1000,F_plot,'LineWidth',3,'Color',[192/255,0,0])
hold on

plot(FdCurve_Sim(:,2)*1000,FdCurve_Sim(:,1),'k','LineWidth',2)
xlim([0 12]);

plot([0,12],[0,0],'k--')

fs = 20;
xlabel('Displacement, d [mm]','FontSize',fs,'FontWeight','bold','FontName','Calibri')
ylabel('Force, F [N]','FontSize',fs,'FontWeight','bold','FontName','Calibri')
set(gca,'FontName','Calibri','FontSize',fs,'FontWeight',...
    'bold');


