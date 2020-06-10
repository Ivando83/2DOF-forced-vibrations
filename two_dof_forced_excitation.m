clear all
close all
clc

%system parameters
%damping of front wheel
c_f = 50;
%stiffness of front wheel
k_f = 70000;
%damping of rear wheel
c_r = 50;
%stiffness of rear wheel
k_r = 70000;
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
p_0 = 0.0;%radians
p_dot_0 = 0;%radians

%excitation force magnitude
force = 2000;%[N]
%frequency of the harmonic force
omega = 1.88*2*pi;%[rad/s]

%eccentricity of the force from the COM of the mass
l_force = 0.5;%[m]

%sampling rate
fs = 100;

%time span
time_span = [0:1/fs:25];

%type of solver
solver = 'anal';

%number of test runs for each computation
n = 10;

switch solver
    case 'dsolve'
        %dsolve
        %symbolic
        for ii = 1:n
            %symbol for independent variable time
            %starting internal loop timer
            tic
            syms t;
            
            %defining symbolic function for verticle bounce
            x = symfun(str2sym('x(t)'),t);
            %velocity
            d1x = diff(x,1);
            %acceleration
            d2x = diff(x,2);
            
            %defining symbolic function for pitch
            p = symfun(str2sym('p(t)'),t);
            %angular velocity
            d1p = diff(p,1);
            d2p = diff(p,2);
            %defining system equations for verticle bounce and pitch
            %verticle bounce
            sysEq1 = (m*d2x + (c_r + c_f)*d1x + (c_f*l_f - c_r*l_r)*d1p + (k_r + k_f)*x + (k_f*l_f - k_r*l_r)*p == force*cos(omega*t));
            %pitch
            sysEq2 = (j*d2p + (c_f*l_f - c_r*l_r)*d1x + (c_f*l_f^2 + c_r*l_r^2)*d1p + (k_f*l_f - k_r*l_r)*x + (k_f*l_f^2 + k_r*l_r^2)*p == l_force*force*cos(omega*t));
            
            
            %solving the two equations
            %the results/output variables are written in alphabetical order
            [p,x] =  dsolve(sysEq1,sysEq2,x(0) == x_0,p(0) == p_0,d1x(0) == x_dot_0,d1p(0) == p_dot_0,t);
            %stopping timer for capturing computational time
            computational_time_sym(ii) = toc;
            %substituting the independent variable t with the input time span
            %bounce
            x_t_sym = double(subs(x,t,time_span));
            %function handle for bounce velocity
            v_fun = matlabFunction(diff(x));
            %bounce velocity
            v_t_sym = feval(v_fun,time_span);
            %pitch motion
            p_t_sym = double(subs(p,t,time_span));
            %function handle for pitch angular velocity
            vp_fun = matlabFunction(diff(p));
            %bounce velocity
            vp_t_sym = feval(vp_fun,time_span);
            %plotting the motions
            figure(1)
            subplot(1,2,1)
            hold on
            aa = plot(time_span,x_t_sym,'-*','color','k')
            bb = plot(time_span,v_t_sym,'-*','color','r')
            legend([aa,bb],'Bounce Motion','Bounce Velocity')
            title('Bounce Motion')
            subplot(1,2,2)
            hold on
            cc = plot(time_span,p_t_sym,'-*','color','k')
            dd = plot(time_span,vp_t_sym,'-*','color','r')
            legend([cc,dd],'Pitch Angular Motion','Pitch Angular Velocity')
            title('Pitch Motion')
        end
        %average computation time for symbolic solution
        sym_time = sum(computational_time_sym)/n;
        
        
    case 'ode45'
        for jj = 1:n
            tic
            %initial conditions vector
            w_0 = [x_0;p_0;x_dot_0;p_dot_0];
            [t,results] = ode45(@(t,w)statefunction(w,t,m,j,k_f,k_r,c_f,c_r,l_f,l_r,force,omega,l_force),time_span,w_0);
            computational_time_num(jj) = toc;
            %bounce motion
            x_t_num = results(:,1);
            %pitch motion
            p_t_num = results(:,2);
            %bounce velocity
            v_t_num = results(:,3);
            %pitch velocity
            vp_t_num = results(:,4);
            
            %plotting the motions
            figure(1)
            subplot(1,2,1)
            hold on
            aa = plot(time_span,x_t_num,'-*','color','k')
            bb = plot(time_span,v_t_num,'-*','color','r')
            legend([aa,bb],'Bounce Motion','Bounce Velocity')
            title('Bounce Motion')
            subplot(1,2,2)
            hold on
            cc = plot(time_span,p_t_num,'-*','color','k')
            dd = plot(time_span,vp_t_num,'-*','color','r')
            legend([cc,dd],'Pitch Angular Motion','Pitch Angular Velocity')
            title('Pitch Motion')
        end
        
        %average computation time for numerical solution
        num_time = sum(computational_time_num)/n;
        
        
    case 'anal'
        for kk = 1:n
            tic
            %system matrices
            %mass matrix
            M = [m,0;0,j];
            %stiffness matrix
            K = [(k_r + k_f),(k_f*l_f - k_r*l_r);(k_f*l_f - k_r*l_r),(k_f*l_f^2 + k_r*l_r^2)];
            %damping matrix
            C = [(c_r + c_f),(c_f*l_f - c_r*l_r);(c_f*l_f - c_r*l_r),(c_f*l_f^2 + c_r*l_r^2)];
            %excitation matrix
            h = [force*cos(omega*time_span);l_force*force*cos(omega*time_span)];%2xlength(time_span) matrix
            
            %for every time step in the time span
            %non homogeneous (particular) solution
            h_star = 0.5*[force;l_force*force];%2x1 matrix
            %frequency response matrix
            %basically a sum of three 2x2 matrices and its inverse
            F_star = inv(-(omega^2)*M + i*omega*C + K);%2x2 matrix
            
            %X star matrix
            X_star = F_star*h_star;%2x1 matrix
            
            %conjugate of X_star
            X_star_bar = conj(X_star);
            
            %particular solution
            %in the following expression
            %2x1 matrix * 1xlength(time_span) matrix
            X_p = X_star*exp(i*omega*time_span) + X_star_bar*exp(-i*omega*time_span);%2xlength(time_span) matrix
            V_p = i*omega*X_star*exp(i*omega*time_span) - i*omega*X_star_bar*exp(-i*omega*time_span);%2xlength(time_span) matrix
            %bounce particular motion
            x_t_p = X_p(1,:);
            v_t_p = V_p(1,:);
            %pitch particular motion
            p_t_p = X_p(2,:);
            vp_t_p = V_p(2,:);
            
            
            %homogeneous solution
            %using polyeigen function to find out the eigenvalues and
            %eigenvectoes of the system
            %in this system, the mass matrix is associated with lambda^2,
            %stiffness matrix is associated with lambda^0 and damping matrix is
            %associated with lambda
            [eigen_vector,lambda] = polyeig(K,C,M);
            
            %the eigen values
            %for this particular case, eigen values are coming out as conjugate
            %complex pairs with negative real parts and thus the system is
            %comparitively weakly damped
            l_1 = lambda(1);
            l_2 = lambda(2);
            l_3 = lambda(3);
            l_4 = lambda(4);
            
            %eigen vector elements
            %first eigen vector
            c_11 = eigen_vector(1,1);
            c_12 = eigen_vector(2,1);
            %second eigenvector
            c_21 = eigen_vector(1,2);
            c_22 = eigen_vector(2,2);
            %3rd eigen vector
            c_31 = eigen_vector(1,3);
            c_32 = eigen_vector(2,3);
            %4th eigenvector
            c_41 = eigen_vector(1,4);
            c_42 = eigen_vector(2,4);
            
            
            
            %modified initial condition matrix
            w_0 = [x_0 - X_star(1) - X_star_bar(1);p_0 - X_star(2) - X_star_bar(2);x_dot_0 - i*omega*X_star(1) + i*omega*X_star_bar(1);p_dot_0 - i*omega*X_star(2) + i*omega*X_star_bar(2)];
            %intermediate eigen vector matrix
            A = [c_11,c_21,c_31,c_41;c_12,c_22,c_32,c_42;l_1*c_11,l_2*c_21,l_3*c_31,l_4*c_41;l_1*c_12,l_2*c_22,l_3*c_32,l_4*c_42];
            %integration constants matrix
            int_const = A\w_0;
            
            %EOM for bounce
            x_t_h = real(c_11*int_const(1)*exp(l_1*time_span) + c_21*int_const(2)*exp(l_2*time_span) + c_31*int_const(3)*exp(l_3*time_span) + c_41*int_const(4)*exp(l_3*time_span));
            v_t_h = real(l_1*c_11*int_const(1)*exp(l_1*time_span) + l_2*c_21*int_const(2)*exp(l_2*time_span) + l_3*c_31*int_const(3)*exp(l_3*time_span) + l_4*c_41*int_const(4)*exp(l_4*time_span));
            %EOM for pitch
            p_t_h = real(c_12*int_const(1)*exp(l_1*time_span) + c_22*int_const(2)*exp(l_2*time_span) + c_32*int_const(3)*exp(l_3*time_span) + c_42*int_const(4)*exp(l_4*time_span));
            vp_t_h = real(l_1*c_12*int_const(1)*exp(l_1*time_span) + l_2*c_22*int_const(2)*exp(l_2*time_span) + l_3*c_32*int_const(3)*exp(l_3*time_span) + l_4*c_42*int_const(4)*exp(l_4*time_span));
            
            
            %bounce motion
            x_t_anal = x_t_h + x_t_p;
            %bounce velocity
            v_t_anal = v_t_h + v_t_p;
            %pitch motion
            p_t_anal = p_t_h + p_t_p;
            %pitch velocity
            vp_t_anal = vp_t_h + vp_t_p;
            
            %computation time for each test run/loop
            computational_time_anal(kk) = toc;
            %plotting
            figure(1)
            subplot(1,2,1)
            hold on
            aa = plot(time_span,x_t_anal,'color','k')
            % bb = plot(time_span,v_t_anal,'-*','color','r')
            % legend([aa,bb],'Bounce Motion','Bounce Velocity')
            title('Bounce Motion')
            subplot(1,2,2)
            hold on
            cc = plot(time_span,p_t_anal,'color','k')
            % dd = plot(time_span,vp_t_anal,'-*','color','r')
            % legend([cc,dd],'Pitch Angular Motion','Pitch Angular Velocity')
            title('Pitch Motion')
            
            
        end
        %average computation time for analytical solution
        anal_time = sum(computational_time_anal)/n;
        
        %only for small or no damping
        eigen_frq_angular = abs(lambda);
        eigen_freq = eigen_frq_angular/(2*pi);
end


