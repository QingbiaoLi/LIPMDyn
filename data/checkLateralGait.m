clear;
clc;
close all;
load LateralGait.txt;

time = LateralGait(:,1);
ycom = LateralGait(:,2);
dycom = LateralGait(:,3);
ddycom = LateralGait(:,4);
ycop = LateralGait(:,5);
PredictTime = LateralGait(:,6);
PredictFP = LateralGait(:,7);
PredictVel= LateralGait(:,8);
dis_ddy = LateralGait(:,9);
fp_wrt_com = LateralGait(:,10);

%%
figure;
f1=figure(1);set(f1,'name','y COM & y COP VS t','numbertitle','off')


% clf; 
hold on; grid on;
plot(time, ycom);
plot(time, ycop,'m');
plot(PredictTime, PredictFP,'go');
xlabel('Time [s]');
ylabel('Position [m]');
legend('y COM','y COP','Foot placement evolve according to time')


%%
figure;
plot(time,dycom);
%%
% figure;
% plot(time,fp_wrt_com);
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
