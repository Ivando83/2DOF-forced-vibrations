clear all
close all
clc

%system parameters
%damping of front wheel
c_f = 100;
%stiffness of front wheel
k_f = 60000;
%damping of rear wheel
c_r = 200;
%stiffness of rear wheel
k_r = 80000;
%mass
m = 1000;
%mass moment of intertia about COM
j = 1000;
%front wheel offset from COM [m]
l_f = 2.5;
%rear wheel offset from COM [m] 
l_r = 2.5;

%initial conditions
%bounce
x_0 = 0.1;
x_dot_0 = 0;
%pitch
p_0 = 0.5;%radians
p_dot_0 = 0;%radians

%base excitation magnitude
base_excitation = 0.1;%[m]

%speed of the car
v = 14;%[m/s], about 50 km/hr

%length of excitation/bump of the road
l_excitation = 5;%[m]

%frequency of the harmonic force
omega = 2*pi*v/l_excitation;%[rad/s]

%phase difference between rear and front excitation
phi = ((l_r + l_f)/l_excitation)*2*pi;


%sampling rate
fs = 100;

%time span
time_span = [0:1/fs:50];

%initial state vector
state_vector_0 = [x_0;p_0;x_dot_0;p_dot_0];

%using ode45 to calculate the state vector

[t, state_vector] = ode45(@(t,state_vector)state_function_base_excitation(t,state_vector,m,j,c_f,c_r,k_f,k_r,l_f,l_r,omega,base_excitation,phi),time_span,state_vector_0);

%bounce motion
x_t = state_vector(:,1);
%pitch motion
p_t = state_vector(:,2);
%bounce velocity
v_t = state_vector(:,3);
%pitch velocity
vp_t = state_vector(:,4);

figure(1)
subplot(1,2,1)
plot(t,x_t)
xlabel('Time[s]')
ylabel('Bounce[m]')

subplot(1,2,2)
plot(t,p_t)
xlabel('Time[s]')
ylabel('Pitch[radians]')


%analytical solution to calculate the amplitude of the particular solution
%(steady state motion)
%frequency response curve to be plotted thus a range of frequencies to be
%taken

%range of speed of car
speed_car = [0:1:200];%km/hr

%range of frequencies
omega_excitation = speed_car*((2*pi)/(l_excitation*3.6));

%mass matrix
M = [m,0;0,j];

%stiffness matrix
K = [k_r + k_f , k_f*l_f - k_r*l_r ; k_f*l_f - k_r*l_r , k_r*(l_r)^2 + k_f*(l_f)^2];

%damping matrix
C = [c_r + c_f , c_f*l_f - c_r*l_r ; c_f*l_f - c_r*l_r , c_r*(l_r)^2 + c_f*(l_f)^2];

for ii = 1:length(omega_excitation)
    
    %the excitation matrix can be written as
    % h = h_c.cos(omega*t) + h_s.sin(omega*t)
    %excitation matrix cosine component
    h_c = [k_r*base_excitation*omega_excitation*sin(phi) + c_f*base_excitation + c_r*base_excitation*cos(phi) ; -k_r*l_r*base_excitation*omega_excitation*sin(phi) + c_f*l_f*base_excitation - c_r*l_r*base_excitation*cos(phi)];
    
    %excitation matrix sine component
    h_s = [-k_f*base_excitation*omega_excitation(ii) - k_r*base_excitation*omega_excitation(ii)*cos(phi) + c_r*base_excitation*sin(phi) ; -k_f*l_f*base_excitation*omega_excitation(ii) + k_r*l_r*base_excitation*omega_excitation(ii)*cos(phi) - c_r*l_r*base_excitation*sin(phi)];
    
    %h star matrix
    h_star = 0.5*(h_c - 1i*h_s);%2x1 matrix
    
    %frequency response matrix
    f_star = inv((-omega_excitation(ii)^2)*M + 1i*omega_excitation(ii)*C + K);%2x2 matrix
    
    X_star = f_star*h_star;%2x1 matrix
    
    x_star = X_star(1);
    p_star = X_star(2);
    
    %cosine component of the bounce motion
    x_c = real(2*x_star);
    %sine component of the bounce motion
    x_s = -imag(2*x_star);
    
    %cosine component of the pitch motion
    p_c = real(2*p_star);
    %sine component of the pitch motion
    p_s = -imag(2*p_star);
    
    %amplitude of the steady state bounce motion
    x_amp(ii) = sqrt((x_c)^2 + (x_s)^2);
    
    %amplitude of the steady state pitch motion
    p_amp(ii) = sqrt((p_c)^2 + (p_s)^2);
end

figure(2)
subplot(1,2,1)
plot(speed_car,x_amp)
xlabel('Car Speed [m/s]')
ylabel('Bounce Amplitude [m]')

subplot(1,2,2)
plot(speed_car,p_amp)
xlabel('Car Speed [m/s]')
ylabel('Pitch Amplitude [radians]')

%speed of car for first resonance mode 
[~,indx1] = max(x_amp);
speed_car_res_mode(1) = speed_car(indx1);
%1st eigen angular frequency (natural angular frequency)
eigen_angular_freq(1) = omega_excitation(indx1);

%speed of car for second resonance mode 
[value,indx2] = max(p_amp);
speed_car_res_mode(2) = speed_car(indx2);
%1st eigen angular frequency (natural angular frequency)
eigen_angular_freq(2) = omega_excitation(indx2);
