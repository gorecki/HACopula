function N = intersectNs(N1, N2)
%INTERSECTNS - Intersects two elements from a nesting semigroup.
%
% N = intersectNs(N1, N2) returns a cell array containing the intersection
% of two elements N1 and N2 from a nesting semigroup. For details, see
% [Gorecki et al., 2017].
%
% References:
% [Gorecki et al., 2017] On Structure, Family and Parameter Estimation
%     of Hierarchical Archimedean copulas. Journal of Statistical Computation 
%     and Simulation, 87(17), 3261-3324
%
%
% Copyright 2018 Jan Gorecki

N = cell(0,0);
for i = 1:length(N1)
    fam1 = N1{i}{1};
    for j = 1:length(N2)
        fam2 = N2{j}{1};
        if strcmp(fam1, fam2)
            r1 = N1{i}{2};
            r2 = N2{j}{2};
            % intersect the intervals
            rmin = max(r1(1), r2(1)); 
            rmax = min(r1(2), r2(2));
            if rmin <= rmax
                % add to the intersection
                N = [N {{fam1, [rmin rmax]}}];
            end
        end
    end
end
        
end
