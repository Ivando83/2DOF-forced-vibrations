clear all %#ok<*CLALL>
close all
clc

%solver
solver = 'anal';

%system parameters
%damping of front wheel
c_f = 100;
%stiffness of front wheel
k_f = 60000;
%damping of rear wheel
c_r = 100;
%stiffness of rear wheel
k_r = 70000;
%mass
m = 1000;
%mass moment of intertia about COM
j = 1000;
%front wheel offset from COM [m]
l_f = 3.5;
%rear wheel offset from COM [m] 
l_r = 2.5;

%initial conditions
%bounce
x_0 = 0.1;
x_dot_0 = 0;
%pitch
p_0 = 1.0;%radians
p_dot_0 = 0;%radians


%position of excitation from the COM acting as the starting point for the
%fmincon solver
length_force_0 = 1.6;

%we are going to find the optimal value of the position of excitation for
%the minimum displacement/vibration in the system
%thus our independent variable in this problem is the position of
%excitation which will be bounded
%lower bound of position of excitation
length_force_lb = -4;%left side of the COM
%upper bound
length_force_ub = 4;%right side of the COM

%amplitude of excitation force
force = 2000;

%angular excitation frequency
omega = 2.4;

%sampling rate
fs = 100;

%time span
time_span = 0:1/fs:75;

%the vibration/aggregate motion is measured with the help of a cost
%function
%the cost function is calculated for each position of excitation force

length_force = length_force_lb:0.01:length_force_ub;

for ii = 1:length(length_force)
    aggregate_motion(ii) = cost_function(length_force(ii),solver,m,j,c_f,c_r,k_f,k_r,l_f,l_r,time_span,force,omega,x_0,p_0,x_dot_0,p_dot_0); %#ok<*SAGROW>
end



%plotting the aggregate motion vs the position of the excitation
figure(1)
plot(length_force,aggregate_motion)
xlabel('Position of Excitation')
ylabel('Aggregate Motion')


%using local search method to find optimum position of excitation for
%minimum excitation/aggregate motion

%for varying angular frequencies
%the local search for optimum value is carried out for each angular
%frequency
for angular_freq = 1:50
    %function handle for cost function
    %below the function handle defines the function for which the
    %independent variable (defined as length_force) will be varied between
    %the defined upper and lower bounds and the amoung the corresponding
    %function values, the minimum function value and the corresponding
    %value of the independent variable will be the output by fmincon
f = @(length_force)cost_function(length_force,solver,m,j,c_f,c_r,k_f,k_r,l_f,l_r,time_span,force,angular_freq,x_0,p_0,x_dot_0,p_dot_0);
length_force_optimum(angular_freq) = fmincon(f,length_force_0,[],[],[],[],length_force_lb,length_force_ub);
end
        

% % %applying global search method for the same
% gs = GlobalSearch;
%
% %defining optimization problem
% problem  =createOptimProblem('fmincon','x0',length_force_0,'objective',f,'lb',length_force_lb,'ub',length_force_ub);
% 
% %running the global search solver
% [length_force_optimum] = run(gs,problem);