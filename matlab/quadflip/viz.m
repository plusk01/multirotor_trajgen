function handle = viz(state, traj, cmd, handle, P)
%VIZ Simulation visualization and drawing tools

if isempty(handle)
    handle{1} = drawEnvironment(state, traj, P);
    handle{2} = drawDesiredAttitude(state, eye(3), []);
    handle{3} = drawDesiredForce(state, zeros(3,1), []);
else
    handle{1} = drawQuad(state, handle{1});
    handle{2} = drawDesiredAttitude(state, cmd.qdes, handle{2});
    handle{3} = drawDesiredForce(state, cmd.Fi, handle{3});
end

end

function hQuad = drawEnvironment(state, traj, P)
    % draw quad in initial state
    hQuad = drawQuad(state, []);
    
    for i = 1:size(traj.p,1)
    
        % only plot at end of segments (dt rate)
        k = mod(i,P.trajDrawPeriod/P.Ts);
        if k ~= 0, continue; end

        % calculate desired attitude from acceleration vector
        qdes = computeAttitude(traj.a(i,:)');
%         qdes = Q.Identity;
        R = Q(qdes).toRotm();

        drawCoordinateAxes(traj.p(i,:)', R, 0.1, 0.2);
    end
    
    % world margin
    s = 1;
    e = size(traj.p,1);
    m = [-1 1 -1 1 -1 1]*0.75; % nominal margin
    datalimits = [min(traj.p(s:e,1)) max(traj.p(s:e,1))...
                  min(traj.p(s:e,2)) max(traj.p(s:e,2))...
                  min(traj.p(s:e,3)) max(traj.p(s:e,3))];
    axis equal; grid on;
    axis(m + datalimits)
    xlabel('X'); ylabel('Y'); zlabel('Z');
    
    % draw coordinate axes
    view(0,0); % force axis to realize its 3d (hack!)
    A = axis;
    k = 0.2;
    O = [A(1);A(3);A(5)];
    drawCoordinateAxes(O, eye(3), k, 1);
end


function handle = drawQuad(state, handle)
    % get rotation of body w.r.t world
    R = Q(state.q).toRotm();
    
    x = state.pos(1);
    y = state.pos(2);
    z = state.pos(3);
    
    x2 = 0.25/2;  % half-distance along x
    y2 = 0.25/2;  % half-distance along y
    z2 = 0.025/2; % half-distance along z
    
    % vertices in body frame
    V = [...
         % top of quad rectangle
         0+x2 0+y2 0+z2;...
         0+x2 0-y2 0+z2;...
         0-x2 0-y2 0+z2;...
         0-x2 0+y2 0+z2;...
         % bottom of quad rectangle
         0+x2 0+y2 0-z2;...
         0+x2 0-y2 0-z2;...
         0-x2 0-y2 0-z2;...
         0-x2 0+y2 0-z2;...
         ];
     
     F = [...
         % top
         1 2 3 4;...
         % side facing +x
         1 5 6 2;...
         % side facing +y
         1 5 8 4;...
         % side facing -x
         3 7 8 4;...
         % side facing -y
         2 6 7 3;...
         % bottom
         5 6 7 8;...
         ];
     
    C = [...
         % top
         0 0 1;...
         % side facing +x
         1 0 0;...
         % side facing +y
         0 1 0;...
         % side facing -x
         0.4 0.4 0.4;...
         % side facing -y
         0.4 0.4 0.4;...
         % bottom
         0.4 0.4 0.4;...
         ];
    
    V = state.pos' + V*R';
     
    if isempty(handle)
        handle = patch('Faces',F,'Vertices',V,...
                       'FaceVertexCData',C,'FaceColor','flat');
    else
        set(handle,'Faces',F,'Vertices',V);
    end
end

function handle = drawDesiredForce(state, F, handle)
    % Scale for aesthetics
    F = 0.05*F;
    
    % TODO: need to rotate?
    
    XX = state.pos(1)+[0 F(1)];
    YY = state.pos(2)+[0 F(2)];
    ZZ = state.pos(3)+[0 F(3)];
    
    if isempty(handle)
        handle = plot3(XX,YY,ZZ);
    else
        set(handle,'XData',XX,'YData',YY,'ZData',ZZ);
    end
end

function handle = drawDesiredAttitude(state, q, handle)
    x = state.pos(1);
    y = state.pos(2);
    z = state.pos(3);

    R = Q(q).toRotm();
    
    handle = drawCoordinateAxes([x;y;z], R', 0.2, 1, handle);
end

function handle = drawCoordinateAxes(O, R, k, alpha, varargin)
%PLOTCOORDINATEFRAME Plot a coordinate frame origin and orientation
%   Coordinate frames have a origin and an orientation. This function draws
%   the coordinate axes in a common frame.

    % k     Size of each axis
    
    if length(varargin)==0 || isempty(varargin{1}), handle = {}; end
    if length(varargin)>=1, handle = varargin{1}; end

    % Create coordinate axes starting at 0
    kk = linspace(0,k,100);
    CX = [kk; zeros(1,length(kk)); zeros(1,length(kk))];
    CY = [zeros(1,length(kk)); kk; zeros(1,length(kk))];
    CZ = [zeros(1,length(kk)); zeros(1,length(kk)); kk];
    
    % First rotate the coordinate frame from the local frame to the
    % orientation of the desired frame in which we want to plot.
    CX = R'*CX;
    CY = R'*CY;
    CZ = R'*CZ;
    
    % Then translate this frame to its origin
    CX = repmat(O, 1, size(CX,2)) + CX;
    CY = repmat(O, 1, size(CY,2)) + CY;
    CZ = repmat(O, 1, size(CZ,2)) + CZ;
    
    % Plot the axes
    ls = '-';
    if isempty(handle)
        handle{1} = plot3(CX(1,:), CX(2,:), CX(3,:),'color','r','linewidth',2,'linestyle',ls);
        handle{2} = plot3(CY(1,:), CY(2,:), CY(3,:),'color','g','linewidth',2,'linestyle',ls);
        handle{3} = plot3(CZ(1,:), CZ(2,:), CZ(3,:),'color','b','linewidth',2,'linestyle',ls);
        
        handle{1}.Color(4) = alpha;
        handle{2}.Color(4) = alpha;
        handle{3}.Color(4) = alpha;
    else
        set(handle{1},'XData',CX(1,:),'YData',CX(2,:),'ZData',CX(3,:));
        set(handle{2},'XData',CY(1,:),'YData',CY(2,:),'ZData',CY(3,:));
        set(handle{3},'XData',CZ(1,:),'YData',CZ(2,:),'ZData',CZ(3,:));
    end
end

function qdes = computeAttitude(a)
% Uses the acceleration vector to reconstruct the desired attitude
% We assume a body flu coordinate frame

    % We add an acceleration in body z to counteract accel due to gravity
    a = a + [0;0;9.81];

    % If the desired acceleration vector is zero, bail
    if sum(a) == 0
        qdes = Q.Identity;
        return;
    end

    % we only care about accel dir
    a = a/norm(a);

    % find the axis of (positive) rotation
    z = [0;0;1];
    axis = cross(a,z);
    axis = axis/norm(axis);
    
    % angle of rotation about the axis
    angle = acos(dot(z,a));
    
    qdes = Q.fromAxisAngle(axis, angle).q;
end

