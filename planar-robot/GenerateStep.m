function [q_d, dq_d] = GenerateStep(t_step, offset, t)
dt = t(2) - t(1);
step = t_step ./ dt;

% First joint
q_d(1,:) = 0*t + offset(1);
q_d(1, step:end) = 0 + offset(1);
dq_d(1,:) = 0*t;

% Second joint
q_d(2,:) = 0*t + offset(2);
q_d(2, step:end) = pi/6 + offset(2);
dq_d(2,:) = 0*t;

end