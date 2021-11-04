function intersections = p_intersectingPoints (prototypes)
% determines the intersection points of the prototypes

% number of prototypes
K = size(prototypes,1);

% there are 4 cases:
%   K <= 1: no intersection points
%   K == 2: infinitely many 'intersection' points
%   K == 3: there are two intersection points on opposite sides
%   K >= 4: normal case

% we can only calculate the intersection points of case K == 3 and the
% normal case
if (K < 3)
    intersections = struct();
    return;
end

% the number of intersections
I = nchoosek (K, 3);
% initialization of interactions cells
intersections(I) = struct();

idx = 1;
% The point where two circles intersect, divides 3 prototypes
for i = 1:K-2
    for j = i+1:K-1
        for k = j+1:K
            % get prototype i, j and k
            w_i = prototypes(i,:);
            w_j = prototypes(j,:);
            w_k = prototypes(k,:);
            
            % get two circles
            circle_ij = circle (w_i, w_j);
            circle_jk = circle (w_j, w_k);
            
            % the intersections points of the two circles is the vector 
            % that is perpendicular to the normal vectors
            cr1 = cross (circle_ij, circle_jk);
            cr1 = cr1/norm(cr1);
            % the second one is on the opposite side
            cr2 = -cr1;
            
            ijk = [i j k];
            % if the 3 closest prototypes of cr are ijk, it's an
            % intersection point we want
            if (check (cr1, prototypes, ijk))
                % store the intersection point
                intersections(idx).point = cr1;
                % and the prototypes it separates
                intersections(idx).prototypes = ijk;
                idx = idx + 1;
            end
            if (check (cr2, prototypes, ijk))
                % store the intersection point
                intersections(idx).point = cr2;
                % and the prototypes it separates
                intersections(idx).prototypes = ijk;
                idx = idx + 1;
            end
        end
    end
end

% remove empty
intersections = intersections(1:idx-1);

end

% Calculate the (normal vector of the) circle that separates two prototypes
function c = circle (w_1, w_2)
    % if the prototypes are on opposite sides, move them a very tiny bit
    if (isequal(w_1+w_2, [0 0 0]))
        w_1 = moveTinyBit (w_1);
        w_2 = moveTinyBit (w_2);
    end

    % calculate two middle points
    mid = (w_1 + w_2)/2;
    % normalize so they are on the sphere
    mid = mid/norm(mid);

    % calculate the vector perpendicular to the prototypes
    p = cross (w_1,w_2);
    % normalize so they are on the sphere
    p = p/norm(p);

    % calculate the normal vector of the circle. Because mid and
    % p both lie on the circle, we can calculate the normal vector 
    % of the circle plane by the cross product of mid and p
    c = cross (p,mid);
    % normalize so they are on the sphere
    c = c/norm(c);
end

% this function checks if the 3 closest prototypes to cr are the same as
% ijk
function same = check (cr, prototypes, ijk)
    % determine the (indices of the) 3 closest prototypes
    [~,so_idx] = sort(vecnorm((cr - prototypes)'));
    so_idx = so_idx(1:3);
    % we need to sort, because ij is also sorted
    so_idx = sort(so_idx);
    same = isequal (so_idx, ijk);
end