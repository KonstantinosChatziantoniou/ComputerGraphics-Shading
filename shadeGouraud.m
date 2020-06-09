function Y = shadeGouraud(p, Vn, Pc, C, S, ka, kd, ks, ncoeff, Ia, I0, X)
% INPUT
%   p   2x3     triangle edge coordinates 2d
%   Vn  3x3     normal vectors of edges 3d
%   Pc  3x1     triangle center of mass
%   ka  3x1     ambient coefs
%   kd  3x1     diffuse coefs
%   ks  3x1     specular coefs
%   ncoeff      scattering level
%   S   3xn     coords of light sources
%   I0  3xn     intensity of each light source
%   X   MxNx3   image with pre-painted triangles
p = p';
clrs = zeros(3,3);
clrs = clrs + ambientLight(ka, Ia);
clrs = clrs + diffuseLight(Pc, Vn, kd, S', I0);
clrs = clrs + specularLight(Pc, Vn, C, ks, ncoeff, S', I0);

[triangle, slopes, clrs, grads] = findActiveEdges(p, clrs');

 % Inital x coord, slope, color and color grad
 x1 = triangle(1,2);
 slope1 = slopes(1);
 c1 = clrs(1,:);
 g1 = grads(1,:);
 x2 = triangle(1,2);
 slope2 = slopes(2);
 c2 = clrs(1,:);
 g2 = grads(2,:);
 for y = (triangle(1,1)):triangle(3,1)
     %Check if we reached the end of the second vertex
     if y == triangle(2,1) %intergrated check: special case -> flat top vertex
         x1 = triangle(2,2);
         slope1 = slopes(3);
         c1 = clrs(2,:);
         g1 = grads(3,:);
         %Check if the third vertex is flat
         if y == triangle(3,1) 
             x2 = triangle(3,2);
             c2 = clrs(3,:);
             % restore x1,c1 if just horizontal line
             if y == triangle(1,1)
                 x1 = triangle(1,2);
                 c1 = clrs(1,:);
             end
         end
     end
     %find the active scan line points and interpolated colors
     sgn = sign(x2-x1);
     if sgn == 0
         sgn = 1;
     end
     points = round(x1):sgn:round(x2);
     c = colorInterp(c1, c2, x1, x2, points);
     X(points,y,1) = c(:,1);
     X(points,y,2) = c(:,2);
     X(points,y,3) = c(:,3);
     %Update the active vertex point
     x1 = x1 + slope1;
     x2 = x2 + slope2;
     %Update the active vertex point color
     c1 = c1 + g1;
     c2 = c2 + g2;
  end
 Y = X;
end

function [tr, slopes, clr, grads] = findActiveEdges(triangle, clr)
    %Sorts the edges based on y axis and calculates slopes.
    % Input :
    %       triangle :      The edges of the triangle
    %       clr :           The color of each edge.
    % Output : 
    %       tr :            Rearranged triangle edges
    %       slopes :        The slopes for each vertex 
    %       clr :           Rearranged edge colors
    %       grads :         The gradient of color on each vertex
    % -----------------------------------------------------------
    triangle = triangle';
    [sortedV, sortedIndex] = sortrows(triangle); %sort triangle points by y
    clr = clr(sortedIndex, :);                     %rearrabnge colors
    grads = zeros(3,3);
    invSlope = [0 0 0];
    pairs = [1 2; 1 3; 2 3];
    for i = 1:3
        grads(i,:) = (clr(pairs(i,1),:) - clr(pairs(i,2),:))/(sortedV(pairs(i,1),1)-sortedV(pairs(i,2),1));
        
        invSlope(i) = (sortedV(pairs(i,1),1)-sortedV(pairs(i,2),1))/(sortedV(pairs(i,1),2)-sortedV(pairs(i,2),2));
        invSlope(i) = 1/invSlope(i);
    end
    tr = sortedV;
    slopes = invSlope;
end

function c = colorInterp(A, B, a, b, x)
    % Finds the color of a point with linear interpolation
    % is VECTORIZED
    % Input:
    %       A :     Color of the first point
    %       a :     x coord of the first point
    %       B :     Color of the last point
    %       b :     x coord of the last point
    %       x :     x coord of the point(s) we want to calculate the color(s)
    % Output:
    %       c :     color(s) of the point(s)
    % -----------------------------------------------------------------------
    a = round(a);
    b = round(b);
    %% check if a b are the same point
    if a == b
        c = A;
        return
    end
    c = A - ((A-B)./(a-b)).*(a - x(:));
end