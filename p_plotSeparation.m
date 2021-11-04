function curves = p_plotSeparation (intersections, prototypes)
% computes the boundary lines between intersection points

% number of prototypes
K = size(prototypes,1);

% there are 4 cases:
%   K <= 1: no boundaries
%   K == 2: one boundary (the great circle between the two prototypes)
%   K == 3: there are three boundaries with the same two intersection
%           points
%   K >= 4: normal case

switch (K)
    case 2
        % get the two prototypes
        w1 = prototypes(1,:);
        w2 = prototypes(2,:);
        
        if (isequal(w1+w2,[0 0 0]))
            w1 = moveTinyBit (w1);
        end
        
        % we need to somehow draw a circle that goes around the whole
        % sphere. There are no intersection points, so we improvise. We
        % calculate three points that are on the circle (no opposites) and
        % at least two on other halves. We use
        % ... the middle point of the two prototypes
        i1 = (w1+w2)/2;
        i1 = i1/norm(i1);
        % ... the point perpendicular to the two prototypes
        i2 = cross(w1,w2);
        i2 = i2/norm(i2);
        % ... and the point in the middle of the previous two calculated
        % point
        i3 = (i1+i2)/2;
        i3 = i3/norm(i3);
        % ... but we can't use this point directly, because then the 3
        % points are on the same half, so we use the opposite point
        i3 = -i3;
        
        intersections(3) = struct();
        % the 3 points should all be connected, so we choose prototypes
        % such that there are two the same between the 3 points. Also, they
        % must be unique, because later we remove duplicates
        intersections(1).point = i1;
        intersections(1).prototypes = [1 2 3];
        intersections(2).point = i2;
        intersections(2).prototypes = [1 2 4];
        intersections(3).point = i3;
        intersections(3).prototypes = [2 3 4];
    case 3
        % get the prototypes
        w1 = prototypes(1,:);
        w2 = prototypes(2,:);
        w3 = prototypes(3,:);
        
        if (isequal(w1+w2, [0 0 0]))
            w1 = moveTinyBit (w1);
            w2 = moveTinyBit (w2);
        end
        if (isequal(w2+w3, [0 0 0]))
            w2 = moveTinyBit (w2);
            w3 = moveTinyBit (w3);
        end
        
        % we can't simply draw two curves between two points
        % that are opposites, because there are infinitely many. We thus
        % use some helpers: the middle points
        i3 = (w1+w2)/2;
        i4 = (w1+w3)/2;
        i5 = (w2+w3)/2;
        
        % make sure the prototypes are such that the three lines can be
        % drawn: p1-i3-p2, p1-i4-p2, p1-i5-p2
        intersections(1).prototypes = [1 2 3 4 5 6];
        intersections(2).prototypes = [7 8 9 10 11 12];
        intersections(3).point = i3;
        intersections(3).prototypes = [1 2 7 8];
        intersections(4).point = i4;
        intersections(4).prototypes = [3 4 9 10];
        intersections(5).point = i5;
        intersections(5).prototypes =  [5 6 11 12];
end

% the number of intersection points
I = length(intersections);
% there may be a boundary between every two intersection points
B = nchoosek(I,2);
boundaries(B) = struct();

idx = 1;
% determine which two intersection points form a boundary
%
% two intersection points separate prototypes m and n if they both separate
% m and n, i.e. if they share m and n
for i = 1:I-1
    for j = i+1:I
        % the prototypes separated by i and j
        prots_i = intersections(i).prototypes;
        prots_j = intersections(j).prototypes;
        % the prototypes shared by the two intersection points
        shared = intersect(prots_i, prots_j);
        if (length(shared) == 2)
            % two prototypes are shared, so there is a boundary
            %
            % store the first intersection point
            boundaries(idx).inter_i = intersections(i).point;
            % ... and the second intersection point
            boundaries(idx).inter_j = intersections(j).point;
            
            idx = idx + 1;
        end
    end
end

% remove empty
boundaries = boundaries(1:idx-1);

% here we calculate the curve points
%
% step is the distance between points on the curve
step = pi/360;
% initialize the curves cells
curves = cell (length(boundaries),1);

% for every boundary, calculate the curve points
for b = 1:length(boundaries)
    inter_i = boundaries(b).inter_i;
    inter_j = boundaries(b).inter_j;
    
    curve = computeCurve (inter_i, inter_j, step);
    % store the curve
    curves(b) = num2cell (curve,[1 2]);
end

end

% compute the curve from inter_i to inter_j (with distances step)
% https://stackoverflow.com/questions/16428444/how-to-find-point-along-arc-in-3d-given-center-start-end-points-radius-ce
function curve = computeCurve (inter_i, inter_j, step)
    % get the normal vector of the great circle P going through intersection
    % point i and j on which the curve lies
    N = cross (inter_i, inter_j);
    % now we need two orthogonal vectors X and Y lying in P
    %
    % inter_i lies in P
    X = inter_i;
    % calculate normalized vector orthogonal to X in P
    Y = cross (N, inter_i)/norm(cross (N, inter_i));
    
    % calculate the angle between the two intersection points
    angle = acos (dot(inter_i,inter_j));
    % ... and how many points there are on the curve
    t = (0:step:angle+step);
    
    % initialize the curve
    curve = zeros(length(t),3);
    % calculate the points on the curve
    for i = 1:length(t)
        curve(i,:) = cos(t(i))*X + sin(t(i))*Y;
    end
end