%% set up geometry
% specify reflection coefficients
beta_wall_freq_domain = -1;
beta_wall_time_domain = -0.9;

% specify constants
c = 1480;
beta_surface = -1;
cutoff_time = 10e-3;
N_max = 200; 
damping = -0.3; % damping coefficient for greens functions

% specify tank geometry

% define tank size
Lx = 0.57;
Ly = 0.34;
Lz = 0.4;

% specify source location
x_source = 0.19;
y_source = 0.17;
z_source = 0.2;
r_source = [x_source; y_source; z_source];

% specify receiver location
x_receiver = 0.36;
y_receiver = 0.17;
z_receiver = 0.11;
r_receiver = [x_receiver; y_receiver; z_receiver];

% create source waveforms
dt = 1e-6; % s
t = (0:dt:20e-3-dt)' - 0.001; %s
t0 = 0; % s
p_source = zeros(length(t),3);

% gaussian pulse
tau = 5e-6;
p_source(:,1) = exp(-((t-t0)./tau).^2);

% decaying sinusoidal pulses
tau = [3 1]*1e-3; % s
f = [1 6]*1e3; % Hz
p_source(:,2:3) = sin(2*pi*f.*(t-t0)).*exp(-(t-t0)./tau);
p_source(t<t0,2:3) = 0;

num_signals = size(p_source,2);

% compute received signal (using fft and greens function)
p_free = compute_time_series_free_field(p_source,dt,r_source,r_receiver,c,0);

% compute frequencies of fft
n = size(p_source,1);
fft_freq = (0:n-1)'./(n*dt);
fft_omega = 2*pi*fft_freq;

% compute fft of source time series
fft_source = fft(p_source);

%% compute tank greens function using images

% g_tank_images = compute_tank_greens_function_images(r_source,r_receiver,fft_omega,Lx,Ly,Lz,c,beta_wall_freq_domain,beta_surface,cutoff_time,damping);

%% compute tank greens function using modes

% g_tank_modes = compute_tank_greens_function_modes(r_source,r_receiver,fft_omega,Lx,Ly,Lz,c,N_max,damping);


%% save or load computationally expensive green's functions
% save('saved_greens_functions.mat','g_tank_modes','g_tank_images');
load('saved_greens_functions.mat','g_tank_modes','g_tank_images');

%% compute received signal from greens functions
tic
fft_receiver_tank_images = g_tank_images.*fft_source;
fft_receiver_tank_modes = g_tank_modes.*fft_source;

p_tank_freq_domain_images = real(ifft(fft_receiver_tank_images));
p_tank_freq_domain_modes = real(ifft(fft_receiver_tank_modes));
toc

%% compute received signal using time domain method
[image_distances,image_coefficients] = compute_source_image_distances_and_reflection_coefficients(r_source,r_receiver,Lx,Ly,Lz,c,beta_wall_time_domain,beta_surface,cutoff_time);
p_tank_time_domain = compute_tank_reflection_time_domain(t,p_source,c,image_distances,image_coefficients);

%% plot results
cmap = flip(cbrewer2('Set1',3));
figure(1);
for i = 1:num_signals
    subplot(num_signals,1,i);
    plot(1e3*t,p_free(:,i),'color','black','linewidth',1,'displayname','Free field');
    hold on;
    plot(1e3*t,p_tank_time_domain(:,i),'color',cmap(1,:),'linewidth',1,'displayname','Time domain images');
    plot(1e3*t,p_tank_freq_domain_images(:,i),'color',cmap(2,:),'linewidth',1,'displayname','Frequency domain images');
    plot(1e3*t,p_tank_freq_domain_modes(:,i),'color',cmap(3,:),'linewidth',1,'linestyle','--','displayname','Frequency domain modes');
    hold off;
    xlim([0 4]);
    xlabel('t (ms)');
    ylabel('p');
end
subplot(3,1,1);
legend('location','northeast');
xlim([0 1]);

set(gcf,'position',[280 120 800 600]);

exportgraphics(gcf,'mode_and_image_full_solution_comparison.png','resolution',600);

% compare greens functions
figure(2);
subplot(2,1,1);
plot(fft_freq,real(g_tank_images));
hold on;
plot(fft_freq,imag(g_tank_images));
hold off;
xlabel('f (Hz)')
ylabel('G');
title('Images');

subplot(2,1,2);
plot(fft_freq,real(g_tank_modes));
hold on;
plot(fft_freq,imag(g_tank_modes));
hold off;
xlabel('f (Hz)');
ylabel('G');
title('Modes');

exportgraphics(gcf,'mode_and_image_greens_function_comparison.png','resolution',600);
