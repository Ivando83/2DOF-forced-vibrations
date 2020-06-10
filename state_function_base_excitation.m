function [dwdt] = state_function_base_excitation(t,w,mass,inertia,damping_front,damping_rear,stiffness_front,stiffness_rear,length_front,length_rear,omega,base_excitation,phi)
w1 = w(1);%bounce motion
w2 = w(2);%pitch motion
w3 = w(3);%bounce velocity
w4 = w(4);%pitch velocity

%bounce velocity
dw1dt = w3;
%pitch velocity
dw2dt = w4;
%bounce acceleration
dw3dt = -((stiffness_front+stiffness_rear)/mass)*w1 - ((damping_front+damping_rear)/mass)*w3 - ((stiffness_front*length_front - stiffness_rear*length_rear)/mass)*w2 - ((damping_front*length_front - damping_rear*length_rear)/mass)*w4 + ((stiffness_front*base_excitation*cos(omega*t))/mass) + ((stiffness_rear*base_excitation*cos(omega*t - phi))/mass) - ((damping_front*omega*base_excitation*sin(omega*t))/mass) - ((damping_rear*omega*base_excitation*sin(omega*t - phi))/mass);

%pitch acceleration
dw4dt = -((stiffness_front*length_front - stiffness_rear*length_rear)/inertia)*w1 - ((damping_front*length_front - damping_rear*length_rear)/inertia)*w3 - ((stiffness_front*(length_front)^2 + stiffness_rear*(length_rear)^2)/inertia)*w2 - ((damping_front*(length_front)^2 + damping_rear*(length_rear)^2)/inertia)*w4 + ((length_front*base_excitation*cos(omega*t)*stiffness_front)/inertia) - ((length_front*damping_front*omega*base_excitation*sin(omega*t))/inertia) - ((length_rear*stiffness_rear*base_excitation*cos(omega*t - phi))/inertia) + ((length_rear*damping_rear*omega*base_excitation*sin(omega*t - phi))/inertia);

%first derivative of the state matrix
dwdt = [dw1dt;dw2dt;dw3dt;dw4dt];
end
