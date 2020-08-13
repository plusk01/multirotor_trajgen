%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Turning Paths Into Trajectories Using Parabolic Blends"
% T. Kunz and M. Stilman, GTech Tech Report, 2011.
%
%
% Parker Lusk
% 12 Aug 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc;

% -------------------------------------------------------------------------
% Path Definition

q = [1 2.5 2 3.5];
vmax = 0.4;
amax = 0.8;
alpha = 0.8; % heuristic "line search" type parameter 0<alpha<1

% -------------------------------------------------------------------------
% Time allocation

if 1
    % Automatic time allocation to turn path into trajectory
    [delT, tb] = solveTiming(q, vmax, amax, alpha);
else
    % Manual specification of trajectory from path
    delT = [1.5 2 1.5]*3;
    tb = [2 4 2 1];
end

% -------------------------------------------------------------------------
% Create Trajectory from Linear Path with Parabolic Blends

n = size(q,2);

% determine the higher-order derivatives at the waypoints
v = diff(q) ./ delT;
a = diff([0 v 0]) ./ tb;

% check blend phase constraint, eq 4
lhs = tb(1:end-1) + tb(2:end);
check = lhs<=2*delT

% actual timing of each linear phase (excludes blend)
tl = delT - lhs./2;

% time of waypoint i
T = zeros(1, n);
for i = 1:n
    T(i) = tb(1)/2 + sum(delT(1:i-1));
end

% total duration of traj
tf = T(n) + tb(n)/2;

% -------------------------------------------------------------------------
% Visualization

% time vector
t = linspace(0,tf,1000);

figure(1), clf;
subplot(311); hold on; grid on;
xlabel('Time [s]'); ylabel('q [m]');
plot(t, qt(t,T,q,v,a,tb));
plot(T, q, 'k'); scatter(T, q, 'k', 'filled');
subplot(312); hold on; grid on;
xlabel('Time [s]'); ylabel('qdot [m/s]');
plot(t, qtdot(t,T,q,v,a,tb));
subplot(313); hold on; grid on;
xlabel('Time [s]'); ylabel('qddot [m/s/s]');
plot(t, qtddot(t,T,q,v,a,tb));

% Error at waypoints:
err = abs(q - qt(T,T,q,v,a,tb))
max(err)

% -------------------------------------------------------------------------
% Helpers

function qq = qt(t, T, q, v, a, tb)

qq = zeros(1,length(t));

lhs1 = T - tb./2;
rhs1 = T + tb./2;
lhs2 = rhs1(1:end-1);
rhs2 = lhs1(2:end);

mask1 = lhs1 <= t' & t' <= rhs1;
mask2 = lhs2 <= t' & t' <= rhs2;
[I1,J1] = find(mask1);
[I2,J2] = find(mask2);

vpad = [0 v];

qq(I1) = q(J1) + vpad(J1).*(t(I1)-T(J1)) + 1/2*a(J1).*(t(I1)-T(J1)+tb(J1)./2).^2;
qq(I2) = q(J2) + v(J2).*(t(I2)-T(J2));


end

function qq = qtdot(t, T, q, v, a, tb)

qq = zeros(1,length(t));

lhs1 = T - tb./2;
rhs1 = T + tb./2;
lhs2 = rhs1(1:end-1);
rhs2 = lhs1(2:end);

mask1 = lhs1 <= t' & t' <= rhs1;
mask2 = lhs2 <= t' & t' <= rhs2;
[I1,J1] = find(mask1);
[I2,J2] = find(mask2);

vpad = [0 v];

qq(I1) = vpad(J1) + a(J1).*(t(I1)-T(J1)+tb(J1)./2);
qq(I2) = v(J2);


end

function qq = qtddot(t, T, q, v, a, tb)

qq = zeros(1,length(t));

lhs1 = T - tb./2;
rhs1 = T + tb./2;
lhs2 = rhs1(1:end-1);
rhs2 = lhs1(2:end);

mask1 = lhs1 <= t' & t' <= rhs1;
mask2 = lhs2 <= t' & t' <= rhs2;
[I1,J1] = find(mask1);
[I2,J2] = find(mask2);

vpad = [0 v];

qq(I1) = a(J1);
qq(I2) = 0;


end

function [delT, tb] = solveTiming(q, vmax, amax, alpha)

m = size(q, 1);
n = size(q, 2);

delT = zeros(1, n-1);
tb = inf(1, n);

for i = 1:length(delT)
    % linear segment timing
    % this linear segment takes as long as the slowest DoF
    delT(i) = maxdof(q(:,i+1), q(:,i), vmax); % eq 21
end

% initialize vel factors
f = ones(1,n-1);

it = 0;

while 1

    % update times of linear segments
    delT = delT./f;

    % create velocities
    v = diff(q) ./ delT;

    vpad = [0 v 0];

    for i = 1:length(tb)
        % blend phase timing
        tb(i) = maxdof(vpad(:,i+1), vpad(:,i), amax);
    end

    % check blend phase constraint, eq 4
    % n.b., we don't care about the rhs of last waypoint, so len(check) == n-1
    lhs = tb(1:end-1) + tb(2:end);
    check = lhs<=2*delT

    if all(check), break; end

    delTpad = [inf delT];
    for i = 1:(n-1)
        if check(i) == 0
            f(i) = alpha*f(i);
%             f(i) = sqrt( min(delTpad(i),delTpad(i+1))./tb(i) );

%             if i>1
%                 f(i-1) = min(f(i-1),f(i));
%             end

            % make sure 0 < f <= 1
            if f(i) <= 0, f(i) = 0.01; end
            if f(i) > 1, f(i) = 1; end
        end
    end

    it = it + 1;
    if it>100, break; end
end
it
end

function T = maxdof(x1, x0, xdot)

n = size(x1, 2);

T = 0;
for j = 1:n
    t = abs(x1 - x0)/xdot;
    if t > T
        T = t;
    end
end
end