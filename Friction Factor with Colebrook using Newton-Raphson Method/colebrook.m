clear;clc;%close all;
%% Main Calculation Section
Re = 1e6; % Reynolds number of the flow

F = @(x,y) (1/sqrt(x)) + 2*log10(y/3.7 + 2.51/(Re*sqrt(x))); % Colebrook equation
dF = @(x,y) (-1/(sqrt(x)*sqrt(x)*sqrt(x)))*(0.5 + 1.09008/(Re*(y/3.7 + 2.51/(Re*sqrt(x))))); % Derivative of Colebrook equation

ed = [0.0,0.0001,0.0005,0.001,0.005,0.01,0.05]; % Relative roughness values
md = [0.011645040997991622,0.013441437692508499,0.01720672984406813,0.019943465840476866,0.030465025820875104,0.037964741876160064,0.0715737538598579]; % Corresponding Friction factor values red from Moody chart

for i = 1:7
    f0 = 1 / ( 3.24*(log(6.9/Re + (ed(i)/3.7)^1.11))^2); % Initial value
    for j = 1:6
        f0 = f0 - F(f0,ed(i))/dF(f0,ed(i)); % Newton-Raphson method
    end
    f0array(i)=f0;
    error(i) = 100*(md(i)-f0array(i))/md(i); % Percentage error calculation
end

%% Plotting Section
figure(1)
plot(ed,md,'o','Color',[0, 0.4470, 0.7410],'LineWidth',2);
hold on;
plot(ed,f0array,'Color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
xlabel('Relative Roughness, \epsilon/D');
ylabel('Friction Factor, f');
legend('Moody Chart','Numerical Colebrook')
grid on;

%% Table Output Section
fg = uifigure;
uit = uitable(fg);
d = {ed(1),f0array(1),md(1),error(1);
    ed(2),f0array(2),md(2),error(2);
    ed(3),f0array(3),md(3),error(3);
    ed(4),f0array(4),md(4),error(4);
    ed(5),f0array(5),md(5),error(5);
    ed(6),f0array(6),md(6),error(6);
    ed(7),f0array(7),md(7),error(7);
};
uit.Data = d;
uit.Position = [20 20 520 400];
uit.ColumnName = {'Relative roughness ($\epsilon/D$)','Friction factor ($f$)','Moody Chart value','Percentage error (%)'};