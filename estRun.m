function [x,y,theta,internalStateOut] = estRun(time, dt, internalStateIn, steeringAngle, pedalSpeed, measurement)
% In this function you implement your estimator. The function arguments
% are:
%  time: current time in [s]
%  dt: current time step [s]
%  internalStateIn: the estimator internal state, definition up to you.
%  steeringAngle: the steering angle of the bike, gamma, [rad]
%  pedalSpeed: the rotational speed of the pedal, omega, [rad/s]
%  measurement: the position measurement valid at the current time step
%
% Note: the measurement is a 2D vector, of x-y position measurement.
%  The measurement sensor may fail to return data, in which case the
%  measurement is given as NaN (not a number).
%
% The function has four outputs:
%  est_x: your current best estimate for the bicycle's x-position
%  est_y: your current best estimate for the bicycle's y-position
%  est_theta: your current best estimate for the bicycle's rotation theta
%  internalState: the estimator's internal state, in a format that can be understood by the next call to this function

% Example code only, you'll want to heavily modify this.
% this needs to correspond to your init function:

% Data from Internal State
B_m = internalStateIn.B;
r_m = internalStateIn.r;
x_m = internalStateIn.x;
y_m = internalStateIn.y;
theta_m = internalStateIn.theta;
P_m = internalStateIn.P;

q_m = [x_m;y_m;theta_m;r_m;B_m];
u_m = [pedalSpeed;steeringAngle];
% V_m = [0.05 0 0.01 0.0001 0;
%        0 0.5 0.0001 0.0001 0;
%        0 0 0.001 0.0001 0.0001;
%        0 0 0 0.001 0;
%        0 0 0 0 0.001];

V_m = [0.05 0 0.01 0.00001 0;
       0 0.3 0.0001 0.00001 0;
       0 0 0.001 0.00001 0.00001;
       0 0 0 0.001 0;
       0 0 0 0 0.001];
   
W_m = [1.0893,1.5333;1.5333,2.9880];

if ~isnan(measurement(1)) & ~isnan(measurement(2))
    % have a valid measurement
    [q_p,A,L] = pred_dyn(q_m,u_m,dt);

    % Prediction
    P_p = A*P_m*A' + L*V_m*L';

    %Measurements
    z = measurement';
    
    % Measurement Update
    [h_m,H,M] = meas_dyn(q_p);
    M = eye(2);
    % Gain
    K = P_p*H'*inv(H*P_p*H'+M*W_m*M');
    q_m = q_p + K*(z-h_m);
    P_m = (eye(5)-K*H)*P_p;
end
% OUTPUTS %%
% Update the internal state (will be passed as an argument to the function
% at next run), must obviously be compatible with the format of
% internalStateIn:

x = q_m(1);
y = q_m(2);
theta = q_m(3);
internalStateOut.x = q_m(1);
internalStateOut.y = q_m(2);
internalStateOut.theta = q_m(3);
internalStateOut.r = q_m(4);
internalStateOut.B = q_m(5);
internalStateOut.P = P_m;
end

%%
function [q_p,A,L] = pred_dyn(q_var,u_var,dt)

% System Model
q_p = [q_var(1)+(5*u_var(1)*q_var(4)*cos(q_var(3)))*dt;
         q_var(2)+(5*q_var(4)*u_var(1)*sin(q_var(3)))*dt;
         q_var(3)+(5*u_var(1)*q_var(4)*tan(u_var(2)))*dt/q_var(5);
         q_var(4);
         q_var(5)];
     
% Estimator
A = [1,0,-5*u_var(1)*q_var(4)*sin(q_var(3))*dt,5*dt*u_var(1)*cos(q_var(3)),0;
     0,1,5*u_var(1)*q_var(4)*cos(q_var(3))*dt,5*dt*u_var(1)*sin(q_var(3)),0;
     0,0,1,(5*dt*u_var(1)*tan(u_var(2)))/q_var(5),-(5*dt*u_var(1)*q_var(4)*tan(u_var(2)))/q_var(5)^2;
     0,0,0,1,0;
     0 0 0 0 1];
 
L = eye(5);

end

function [h_m,H,M] = meas_dyn(q_var)

% Measurement Model
h_m = [q_var(1) + q_var(5)*cos(q_var(3))/2;q_var(2) + q_var(5)*sin(q_var(3))/2];
     
% Estimator

H = [1,0,-(q_var(5)*sin(q_var(3)))/2,0,cos(q_var(3))/2;
     0,1,(q_var(5)*cos(q_var(3)))/2,0,sin(q_var(3))/2];

 M = eye(2);
end