%plotting the frequency repsonse graph
clear all %#ok<CLALL>
close all
clc

%state parameters
%damping of front wheel
c_f = 100;
%stiffness of front wheel
k_f = 60000;
%damping of rear wheel
c_r = 100;
%stiffness of rear wheel
k_r = 60000;
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
p_0 = 1.0;%radians
p_dot_0 = 0;%radians

%excitation force magnitude
force = 2000;%[N]
%frequency range of the harmonic force
omega = 0:0.1:50;%[rad/s]

%eccentricity of the force from the COM of the mass
l_force = 1.6;%[m]

%sampling rate
fs = 100;

%time span
time_span = 0:1/fs:75;

%type of solver
solver = 'anal';


switch solver
    case 'anal'
        for kk = 1:length(omega)
            tic
        %system matrices
        %mass matrix
        M = [m,0;0,j];
        %stiffness matrix
        K = [(k_r + k_f),(k_f*l_f - k_r*l_r);(k_f*l_f - k_r*l_r),(k_f*l_f^2 + k_r*l_r^2)];
        %damping matrix
        C = [(c_r + c_f),(c_f*l_f - c_r*l_r);(c_f*l_f - c_r*l_r),(c_f*l_f^2 + c_r*l_r^2)];
        %excitation matrix
        h = [force*cos(omega(kk)*time_span);l_force*force*cos(omega(kk)*time_span)];%2xlength(time_span) matrix 
        
        %for every time step in the time span
        %non homogeneous (particular) solution
        h_star = 0.5*[force;l_force*force];%2x1 matrix
        %frequency response matrix
        %basically a sum of three 2x2 matrices and its inverse
        F_star = inv(-(omega(kk)^2)*M + 1i*omega(kk)*C + K);%2x2 matrix
        
        %X star matrix
        X_star = F_star*h_star;%#ok<*MINV> %2x1 matrix
        
        %conjugate of X_star
        X_star_bar = conj(X_star);
        
        %cos and sin coefficients of particular solution
        %bounce
        x_cos = real(2*X_star(1));
        x_sin = -imag(2*X_star(1));
        
        %pitch
        p_cos = real(2*X_star(2));
        p_sin = -imag(2*X_star(2));
        
        %amplitude of the motion for the particular omega
        %bounce
        x_amp(kk) = sqrt((x_cos)^2 + (x_sin)^2);
        %pitch
        p_amp(kk) = sqrt((p_cos)^2 + (p_sin)^2);
        
        

        end

figure(1)
subplot(1,2,1)
plot(omega,x_amp)
title('Bounce Frequency Response')
xlabel('Angular Frequency')
ylabel('Bounce Amplitude')

subplot(1,2,2)
plot(omega,p_amp)
title('Pitch Frequency Response')
xlabel('Angular Frequency')
ylabel('Pitch Amplitude')

%static vertical deflection incase of a static load
x_st = x_amp(1);

%static angular displacment incase of a static load
p_st = p_amp(1);


end