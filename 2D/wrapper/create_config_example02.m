%% WRAPPER FOR FASTCYCLES 
clear all
close all


% Chose the number of points
nb_element = 2^12;

% Create structure to modify
sim = setdefaults_v01(nb_element);


%% Hyperparameters
% Name of the simulation
sim.simulation_name = 'example02';

% Number of openmp threads
sim.omp_threads = 4     ;  

% Control of the frequency dependent cut
sim.rho_c = 200.0; % No need to change this parameter usually

% Control the first term Tw(1) = eta_w * L_fault / cs
sim.eta_w = 1; % No need to change this parameter usually

% Control on time discretisation 
sim.beta_min = 0.25; 
% If the simulation is not converging, or the 
% time accuracy is not good enough (for fully-dynamic simulation), this 
% parameter can be decreased

% Control on Newton Raphson 
sim.nr_accuracy = 1e-6;

% Control on the time solver
sim.tol_solver = 1e-4;

% Quasi-dynamic of fully dynamic
sim.quasi_dynamic =1;

%% Output 
sim.stride_time = 1; % Only save data every 'sim.stride_time' time step

sim.freq_writing_file = 1000; 
% Write files every sim.freq_writing_file x sim.stride_time
% the above value has to be smaller than sim.max_it.
% sim.freq_writing_file has also a strong influence of the memory usage of
% the program, it can be reduced to reduce memory usage.

sim.max_it = 7000; % the maximum number of iteration


%% Variables
% Main
% Friction
sim.friction_law = 'RateState_R'; % Only this version is available
sim.evolution_law = 'aging'; % Only this version is available

% Material and loading
sim.mu = 40e9; % Shear modulus
sim.cs = 3.464e3; % S-wave velocity
sim.cp = 6.e3; % P-wave velocity

%% Create mask
sim.mask(:) = 1;
sim.mask(1:2^10+2^9) = 0;
sim.mask(end-2^10-2^9+1:end) = 0;

%% Friction
sim.a(:) = 0.012;
sim.b(:) = 0.015;
sim.Dc(:) = 0.01;
sim.f0(:) = 0.6;
sim.V0(:) = 1e-9;

%% Geometry
% Total length of the fault
sim.L_fault = 4*0.1024E+05; 

% Create x and ds
ds = sim.L_fault/sim.nb_element;
x = -sim.L_fault/2+ds/2:ds:sim.L_fault/2;

% Amplitude to wavelength ratio
alpha = 1e-4;

% Parameter for generation of rough geometry
aa=10000000;         % [m], big value make the Von Karman PSD tends to the one in Dunham et al., 2011b (autocorrelation distance)
epsi= 2*sqrt(2)*pi*aa*alpha;   % Theoretical value of epsilon based on alpha and aa to match Dunham et al., 2011b

% Order of the von Karman
kap=1; % For this value, the Von Karman PSDF asymtotically converge toward P(k)=k^(-3) which is self-similar   

% Minimum wavelength of the roughness
lambda_min = 200;

% Random number for the phase
seedXi=squeeze(rand(1,sim.nb_element));

% Generate geometry
sim.height = von_Karman(epsi,aa,kap,seedXi,ds,lambda_min);

figure(1)
plot(x/1000,sim.height)
xlabel('distance along fault (km)')



%% Initial parameters
sim.sigmaN(:) = -100e6;
sim.V(:) = 1e-9;
sim.theta(:) = 0.0;

%% Loading
sim.shear_loading_dot(:) = 0.01; % Add shear traction rate loading on the fault
sim.normal_loading_dot(:) = 0.0; % Add normal traction rate loading on the fault
sim.V_mask = 0.0; % block the fault, where sim.mask = 0 

%% Write config file
directory_mkdir = ['../problems/' sim.simulation_name];
mkdir(directory_mkdir)
writeconfigfile_v01(['../problems/' sim.simulation_name '/'],sim)


%%
% Compute Lb
Lb = sim.mu*2*(1-(sim.cs/sim.cp)^2)*sim.Dc(1)/(abs(sim.sigmaN(1))*sim.b(1));
            
% Compute Lnuc from Viesca [2016]
Lnuc = 2*Lb / (pi*(1-sim.a(1)/sim.b(1))^2);




%% Function that generate fractal (from Tatsuhiko Saito)
function [Xi]=von_Karman(epsi,ax,kap,seedXi,dx,lambda_min)
%
%           ver 1.0   2004.12.27
%           T. Saito, modified by P. Romanet
%                                                 
%

% Lenght of the array
nx=length(seedXi);

% Create wavenumber
kx = 2*pi/(nx*dx)*[0:nx/2-1, -nx/2:-1];

% Uniformely distributed phase between [0:2pi]
pang=2*pi*seedXi;

% Von Karman (Sato et al., 2012, p 24)
yy_c=2*pi^(1/2)*gamma(kap+1/2)*epsi^2*ax./(gamma(kap).*(1+(kx.*ax).^2).^(kap+0.5));
% yy_c1= (2*pi)^3*(1e-4)^2./kx.^3;

% Realisation of random media (Sato et al., 2012, p 21)
Xi_kx=sqrt(yy_c).*exp(1i.*pang)*sqrt(nx*dx);      % [m*m]

% Remove mean of the signal
Xi_kx(1) = 0.0;

% Remove small perturvation
% Calculate max wavenumber
kmax = 2*pi/lambda_min;

% Remove wavenumber that are higher than max wavenumber
index = find(abs(kx)>kmax);
Xi_kx(index) = 0;

% ---   Spectrum to Fluctuation
Xi=ifft(Xi_kx)/dx;  
Xi = real(Xi);


% Calculate RMS
alpha_check = sqrt(dx*sum(Xi.^2)/((nx*dx)))/(nx*dx);
ratio = epsi/(2*sqrt(2)*pi*ax)/alpha_check;
disp(['Theoretical alpha =',num2str(epsi/(2*sqrt(2)*pi*ax)),' Calculated alpha =',num2str(alpha_check)])
disp(['Ratio=',num2str(ratio)])

end





