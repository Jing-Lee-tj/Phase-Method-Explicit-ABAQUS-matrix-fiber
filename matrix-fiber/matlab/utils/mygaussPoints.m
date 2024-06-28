function [points, weights] = mygaussPoints(n,ndim)
    if ndim==1
        [points, weights] = gaussPoints1D(n);
    elseif ndim==2
        [points, weights] = gaussPoints2D(n);
    elseif ndim==3
        [points, weights] = gaussPoints3D(n);
    end
end

function [points, weights] = gaussPoints1D(n)
    % This function computes the Gauss-Legendre points and weights for 1D integration
    % n: number of points
    % points: Gauss points in the interval [-1, 1]
    % weights: corresponding weights
    
    % Legendre polynomial of degree n
    P = legendreP(n, sym('x'));
    
    % Roots of the Legendre polynomial are the Gauss points
    points = double(vpasolve(P == 0));
    
    % Weights calculation
    weights = zeros(n, 1);
    for i = 1:n
        weights(i) = 2 / ((1 - points(i)^2) * (subs(diff(P), points(i)))^2);
    end
end

function [points, weights] = gaussPoints2D(n)
    % This function computes the Gauss-Legendre points and weights for 2D integration
    % n: number of points in each dimension
    % points: Gauss points in the interval [-1, 1]x[-1, 1]
    % weights: corresponding weights
    
    [points1D, weights1D] = gaussPoints1D(n);
    
    % Initialize points and weights
    points = zeros(n^2, 2);
    weights = zeros(n^2, 1);
    
    % Compute 2D points and weights
    index = 1;
    for i = 1:n
        for j = 1:n
            points(index, :) = [points1D(i), points1D(j)];
            weights(index) = weights1D(i) * weights1D(j);
            index = index + 1;
        end
    end
end

function [points, weights] = gaussPoints3D(n)
    % This function computes the Gauss-Legendre points and weights for 3D integration
    % n: number of points in each dimension
    % points: Gauss points in the interval [-1, 1]x[-1, 1]x[-1, 1]
    % weights: corresponding weights
    
    [points1D, weights1D] = gaussPoints1D(n);
    
    % Initialize points and weights
    points = zeros(n^3, 3);
    weights = zeros(n^3, 1);
    
    % Compute 3D points and weights
    index = 1;
    for i = 1:n
        for j = 1:n
            for k = 1:n
                points(index, :) = [points1D(i), points1D(j), points1D(k)];
                weights(index) = weights1D(i) * weights1D(j) * weights1D(k);
                index = index + 1;
            end
        end
    end
end
