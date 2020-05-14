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
x_m = internalStateIn.x;
y_m = internalStateIn.y;
theta_m = internalStateIn.theta;
r_m = internalStateIn.r;
B_m = internalStateIn.B;
P_m = internalStateIn.P;

q_m = [x_m;y_m;theta_m;r_m;B_m];
u_m = [pedalSpeed;steeringAngle];
V_m = [0.6 0 0.1 0.001 0;
       0 0.6 0.1 0.001 0;
       0 0 0.6 0.001 0.001;
       0 0 0 0.05 0;
       0 0 0 0 0.05];
W_m = [1.0893,1.5333;1.5333,2.9880];


if ~isnan(measurement(1)) & ~isnan(measurement(2))
    % have a valid measurement
    % Sigma Point
    sigma = zeros(5,10);
    sigma(:,1:5) = q_m + sqrt(5*P_m);
    sigma(:,6:10) = q_m - sqrt(5*P_m);

    % Prediction Sigma Points
    qp_sigma = zeros(5,10);
    q_p = 0;
    for i = 1:10
        qp_sigma(:,i) = pred_dyn(sigma(:,i),u_m,dt);
        q_p = q_p + qp_sigma(:,i);
    end
    q_p = (1/10)*q_p;
    P_p = 0;
    for i=1:10
        P_p = P_p + (qp_sigma(:,i)-q_p)*(qp_sigma(:,i)-q_p)';
    end
    P_p = real((1/10)*P_p + V_m);

    % Measurement Update
    z_sigma = zeros(2,10);
    z_hat = 0;
    for i = 1:10
        z_sigma(:,i) = meas_dyn(qp_sigma(:,i));
        z_hat = z_hat + z_sigma(:,i);
    end
    z_hat = (1/10)*z_hat;

    P_zz = 0;
    P_qz = 0;
    for i=1:10
        P_zz = P_zz + (z_sigma(:,i)-z_hat)*(z_sigma(:,i)-z_hat)';
        P_qz = P_qz + (qp_sigma(:,i)-q_p)*(z_sigma(:,i)-z_hat)';
    end
    P_zz = real((1/10)*P_zz + W_m);
    P_qz = real((1/10)*P_qz);

    z = measurement';

    K = P_qz*inv(P_zz);
    q_m = real(q_p + K*(z-z_hat));
    P_m = real(P_p - K*P_zz*K');
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
function q_p = pred_dyn(q_var,u_var,dt)

% System Model
q_p = [q_var(1)+(5*u_var(1)*q_var(4)*cos(q_var(3)))*dt;
         q_var(2)+(5*q_var(4)*u_var(1)*sin(q_var(3)))*dt;
         q_var(3)+(5*u_var(1)*q_var(4)*tan(u_var(2)))*dt/q_var(5);
         q_var(4);
         q_var(5)];
     
end

function h_m = meas_dyn(q_var)

% Measurement Model
h_m = [q_var(1) + q_var(5)*cos(q_var(3))/2;q_var(2) + q_var(5)*sin(q_var(3))/2];

end