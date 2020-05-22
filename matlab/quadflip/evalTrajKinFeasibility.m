function feasible = evalTrajKinFeasibility(traj, t, P)
%EVALTRAJDYNFEASIBILITY Determine if each segment is dynamically feasible
%   Detailed explanation goes here

% number of trajectory segments
Nsegs = length(t)-1;

% keep track of the feasibility of each segment
feasible = zeros(1,Nsegs);

% keep track of segment start index
start = 1;

for s = 1:Nsegs
    
    seg.v = traj.v(start:traj.sidx(s), :);
    seg.a = traj.a(start:traj.sidx(s), :);
    seg.j = traj.j(start:traj.sidx(s), :);
    seg.s = traj.s(start:traj.sidx(s), :);
    
    feasible(s) = all(all(abs(seg.v) < P.vmax)) & ...
                    all(all(abs(seg.a) < P.amax)) & ...
                    all(all(abs(seg.j) < P.jmax)) & ...
                    all(all(abs(seg.s) < P.smax));
    
    % cue up next start index
    start = traj.sidx(s) + 1;
end


end