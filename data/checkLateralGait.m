clear;
clc;
close all;
load LateralGait.txt;

time = LateralGait(:,1);
xcom = LateralGait(:,2);
dxcom = LateralGait(:,3);
ddxcom = LateralGait(:,4);
xcop = LateralGait(:,5);

xcom_e_f = LateralGait(:,6);
dxcom_e_f = LateralGait(:,7);
xcop_e_f= LateralGait(:,8);
xcapture2_e_f_predict = LateralGait(:,9);
xcapture2_e_f_observe = LateralGait(:,10);

fp_wrt_com = LateralGait(:,11);
com_wrt_stance_foot = LateralGait(:,12);

coeff_v0 = LateralGait(:,13);
coeff_vd = LateralGait(:,14);
coeff_com_offset= LateralGait(:,15);

%%
%{
result(1) = figure('name','online gain adpation for bipedal walking','units','normalized','outerposition',[0 0 1 1]);

% subplot(2,1,1)
plot(time,coeff_v0,'g--','LineWidth',3);
hold on; grid on;
title('online gain adpation for bipedal walking','FontSize',30);
plot(time,coeff_vd,'b-','LineWidth',3);
 plot(time, coeff_com_offset,'r-.','LineWidth',3);

xlabel('Time (s)','FontSize',30);
ylabel('Coefficients (Sagittal)','FontSize',30);
% legend('Offset coefficient','CoM Pos coefficient',' CoM Vel coefficient','location','southeast')
set(gca,'fontsize',30,'FontWeight','bold');
% legend({'Current Velocity coeff','Desired Velocity coeff'},'FontSize',30,'Location','northwest')
legend({'Desired Velocity coeff','Current Velocity coeff','CoM offset'},'FontSize',30,'Location','northwest')

% saveimages_name =sprintf('%s/online_gain_adpation_for_bipedal_walking',saveimages_dir);
% saveas(result(1),saveimages_name);  
%  print(result(1),saveimages_name,'-dpng','-r300');
% print(result(1),saveimages_name,'-dpdf','-r0');
%}
%% plot result -  CoM vel and desired vel(Sagittal)
result(2)=figure('name','Validation of online estimation method (Sagittal vel)','units','normalized','outerposition',[0 0 1 1]);

plot(time,dxcom,'k','LineWidth',3);
hold on; grid on;
plot(time,dxcom_e_f,'r','LineWidth',3);
title('Sagittal data in support foot frame','FontSize',25);
target_vel= ones(size(dxcom,1),1)*0.3;
plot(time,target_vel,'g','LineWidth',3);
% Reference of Global COM pos
% plot(time, com_ref_x,'g--');

xlabel('Time (s)','FontSize',25);
ylabel('CoM Velocity (sagittal) (m/s)','FontSize',30);
legend({'Global CoM vel',' Target Vel'},'FontSize',30,'Location','southeast')
set(gca,'fontsize',30,'FontWeight','bold');
% saveas(result(3),'result/images/Validation of online estimation method (Local Sagittal).png');  
% saveimages_name =sprintf('%s/Validation of online estimation method (Vel Sagittal)',saveimages_dir);
% saveas(result(3),saveimages_name); 
% print(result(2),saveimages_name,'-dpng','-r300');

%%
%{
figure;
plot(time,fp_wrt_com);
%%
figure;
plot(time,com_wrt_stance_foot);
result(3)=figure('name','Validation of online estimation method (Sagittal vel)','units','normalized','outerposition',[0 0 1 1]);

plot(time,xcom,'k','LineWidth',3);
hold on; grid on;
plot(time,xcom_e_f,'r','LineWidth',3);
title('Sagittal data in support foot frame','FontSize',25);
target_vel= ones(size(dxcom,1),1)*0.5;
% plot(time,target_vel,'g','LineWidth',3);
% Reference of Global COM pos
% plot(time, com_ref_x,'g--');

xlabel('Time (s)','FontSize',25);
ylabel('CoM Velocity (sagittal) (m/s)','FontSize',30);
legend({'Global CoM vel',' Target Vel'},'FontSize',30,'Location','southeast')
set(gca,'fontsize',30,'FontWeight','bold');
% saveas(result(3),'result/images/Validation of online estimation method (Local Sagittal).png');  
% saveimages_name =sprintf('%s/Validation of online estimation method (Vel Sagittal)',saveimages_dir);
% saveas(result(3),saveimages_name); 
% print(result(2),saveimages_name,'-dpng','-r300');
%}
%{
figure;
% return;
% f2=figure(2);set(f2,'name','PredictFP VS PredictTime + dyCOM vs t','numbertitle','off')
subplot(2,2,1);
% clf;
hold on; grid on;
plot(time, dycom);
plot(PredictTime, PredictVel,'ro');
xlabel('PredictTime [s]');
ylabel('PredictFP [m]');
% set(f2,'name','PredictFP VS PredictTime + dyCOM vs t','numbertitle','off')



% f3=figure(3);set(f3,'name','PredictFP  VS PredictTime','numbertitle','off')

subplot(2,2,2);
% clf;
hold on; grid on;
% plot(time, dycom);
plot(PredictTime, PredictFP);
% plot(time, PredictVel,'ro');
plot(PredictTime, PredictFP,'ro');
xlabel('PredictTime [s]');
ylabel('PredictFP [m]');

% return;
% f4=figure(4); set(f4,'name','dyCOM VS t','numbertitle','off')
subplot(2,2,3);
% clf;
hold on; grid on;
plot(time, dycom);
xlabel('Time [s]');
ylabel('dycom [m]');


% f5=figure(5); 
% figure;
% clf;
subplot(2,2,4);
hold on; grid on;
plot(time, dis_ddy);
% set(f5,'name','dis_ddy VS t','numbertitle','off')
title('dis_ddy VS t')

%}
