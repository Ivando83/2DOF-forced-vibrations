function [dwdt] = statefunction(w,t,mass,mass_inertia,stiffness_front,stiffness_rear,damping_front,damping_rear,offset_front,offset_rear,force_excitation,omega,offset_excitation)

w1 = w(1);%bounce motion
w2 = w(2);%pitch motion
w3 = w(3);%bounce velocity
w4 = w(4);%pitch velocity

dw1dt = w3;%bounce velocity
dw2dt = w4;%pitch velocity
%bounce acceleration
dw3dt = (force_excitation*cos(omega*t)/mass)-((damping_rear + damping_front)/mass)*w3 - ((damping_front*offset_front - damping_rear*offset_rear)/mass)*w4 - ((stiffness_rear + stiffness_front)/mass)*w1 - ((stiffness_front*offset_front - stiffness_rear*offset_rear)/mass)*w2;
%pitch acceleration
dw4dt = (offset_excitation*force_excitation*cos(omega*t)/mass_inertia)-((damping_front*offset_front - damping_rear*offset_rear)/mass_inertia)*w3 - ((damping_rear*(offset_rear^2) + damping_front*(offset_front)^2)/mass_inertia)*w4 - ((stiffness_front*(offset_front)^2 + stiffness_rear*(offset_rear)^2)/mass_inertia)*w2 - ((stiffness_front*offset_front - stiffness_rear*offset_rear)/mass_inertia)*w1;

%first order derivative of the state vector
dwdt = [dw1dt;dw2dt;dw3dt;dw4dt];
end