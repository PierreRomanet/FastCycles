%% Read data
clear all
close all
clc

% Chose the directory of the data to read
direc = '../problems/example02/';

% Import the velocity, theta and shear traction
f = loadandprocessdata(direc,'V','theta','St','Sn');
input = load_input(direc);


% Create vector of distance 
L = 4*0.1024E+05;
ds = L/size(f.V,2);
x = -L/2+ds/2:ds:L/2;


%% Velocity vs time step

% Create figure
figure('position',[1 1 1000 1000])

% Plot the result
pcolor(x/1000,1:size(f.V,1),log10(f.V));

% Decoration
% Labels
xlabel('Position along the fault (km)')
ylabel('Time step')

% Colorbar
clb = colorbar;
caxis([-11 10])
clb.Label.String = 'Slip velocity (m/s)';
for j=1:length(clb.TickLabels)
    clb.TickLabels{j} = ['10^{' clb.TickLabels{j} '}'];
end

% General properties
set(gca,'Fontsize',22)
set(gcf,'color','w')


shading interp


%% Elastic Stress evolution along the fault
figure('position',[1 1 1000 1000])
for i=1:10:length(f.time)
    clf
    subplot(2,1,1)
    % Plot
    plot(x,f.St(i,:)/1e6)

    % Title
    title(['Velocity at time: ' seconds2duration(f.time(i))])

    % Limit label
    xlim([min(x) max(x)])
    ylim([min(f.St,[],"all")*1.1/1e6 max(f.St,[],"all")*1.1/1e6])

    % Labels
    xlabel('Position along the fault (km)')
    ylabel('Elastic shear traction (MPa)')

    % subplot properties
    set(gca,'Fontsize',20)

    subplot(2,1,2)
    % Plot
    plot(x,f.Sn(i,:)/1e6)

    % Limit label
    xlim([min(x) max(x)])
    ylim([min(f.Sn,[],"all")*1.1/1e6 max(f.Sn,[],"all")*1.1/1e6])

    % Labels
    xlabel('Position along the fault (km)')
    ylabel('Elastic normal traction (MPa)')

    % subplot properties
    set(gca,'Fontsize',20)

    % General properties
    set(gcf,'color','w')

    % Pause
    pause(0.1)
end


%% Total shear traction
%% Elastic Stress evolution along the fault
figure('position',[1 1 1000 1000])
for i=1:10:length(f.time)
    clf

    % Calculate total shear traction
    St_tot = -input.Sn_ini.*input.a.*asinh(squeeze(f.V(i,:))./(2.*input.V0).*exp(input.f0./input.a+input.b./input.a.*squeeze(f.theta(i,:))));
    St_tot = St_tot(input.mask==1);

    % Plot
    plot(x(input.mask==1),St_tot/1e6)

    % Title
    title(['Velocity at time: ' seconds2duration(f.time(i))])

    % Limit label
    xlim([min(x(input.mask==1)) max(x(input.mask==1))])
%     ylim([min(St_tot,[],"all")*1.1/1e6 max(St_tot,[],"all")*1.1/1e6])

    % Labels
    xlabel('Position along the fault (km)')
    ylabel('Total shear traction (MPa)')

    % subplot properties
    set(gca,'Fontsize',20)

    % General properties
    set(gcf,'color','w')

    % Pause
    pause(0.1)
end
