% Konstantinos Chatziantoniou 8941 konstantic@ece.auth.gr
% Aristotle University of Thessaloniki
% Computer Graphics
% 3rd Assignment - 2020/06/19
function Y = shadePhong(p, Vn, Pc, C, S, ka, kd, ks, ncoeff, Ia, I0, X)
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
   
    
    [triangle, slopes, ka, ka_grads, kd, kd_grads, ks, ks_grads, Vn, Vn_grads] = findActiveEdges(p', ka', kd', ks', Vn');
    
     % Inital x coord, slope, color and color grad
     x1 = triangle(1,2);
     slope1 = slopes(1);
     %c1 = clrs(1,:);
     ka1 = ka(1,:); kd1 = kd(1,:); ks1 = ks(1,:);
     vn1 = Vn(1,:); 
     %g1 = grads(1,:);
     gka1 = ka_grads(1,:);gkd1 = kd_grads(1,:);gks1 = ks_grads(1,:);
     gvn1 = Vn_grads(1,:);
     x2 = triangle(1,2);
     slope2 = slopes(2);
     %c2 = clrs(1,:);
     ka2 = ka(1,:); kd2 = kd(1,:); ks2 = ks(1,:);
     vn2 = Vn(1,:); 
     %g2 = grads(2,:);
     gka2 = ka_grads(2,:);gkd2 = kd_grads(2,:);gks2 = ks_grads(2,:);
     gvn2 = Vn_grads(2,:);
     for y = (triangle(1,1)):triangle(3,1)
         %Check if we reached the end of the second vertex
         if y == triangle(2,1) %intergrated check: special case -> flat top vertex
            x1 = triangle(2,2);
            slope1 = slopes(3);
            %c1 = clrs(2,:);
            ka1 = ka(2,:); kd1 = kd(2,:); ks1 = ks(2,:);
            vn1 = Vn(2,:); 
            %g1 = grads(3,:);
            gka1 = ka_grads(3,:);gkd1 = kd_grads(3,:);gks1 = ks_grads(3,:);
            gvn1 = Vn_grads(3,:);
            %Check if the third vertex is flat
            if y == triangle(3,1) 
                x2 = triangle(3,2);
                %c2 = clrs(3,:);
                ka2 = ka(3,:); kd2 = kd(3,:); ks2 = ks(3,:);
                vn2 = Vn(3,:); 
                % restore x1,c1 if just horizontal line
                if y == triangle(1,1)
                    x1 = triangle(1,2);
                    %c1 = clrs(1,:);
                    ka1 = ka(1,:); kd1 = kd(1,:); ks1 = ks(1,:);
                    vn1 = Vn(1,:); 
                end
            end
         end
         %find the active scan line points and interpolated colors
         sgn = sign(x2-x1);
         if sgn == 0
             sgn = 1;
         end
         points = round(x1):sgn:round(x2);
         c = colorInterp(Pc, S,C, ncoeff, Ia, I0, ka1, ka2, kd1, kd2, ks1, ks2, vn1, vn2, x1, x2, points);
         X(points,y,1) = c(:,1);
         X(points,y,2) = c(:,2);
         X(points,y,3) = c(:,3);
         %Update the active vertex point
         x1 = x1 + slope1;
         x2 = x2 + slope2;
         %Update the active vertex point color
         %c1 = c1 + g1;
         ka1 = ka1 + gka1;
         kd1 = kd1 + gkd1;
         ks1 = ks1 + gks1;
         vn1 = vn1 + gvn1;
         %c2 = c2 + g2;
         ka2 = ka2 + gka2;
         kd2 = kd2 + gkd2;
         ks2 = ks2 + gks2;
         vn2 = vn2 + gvn2;
      end
     Y = X;
    end
    
    function [tr, slopes, ka, ka_grads, kd, kd_grads, ks, ks_grads, Vn, Vn_grads] = findActiveEdges(triangle, ka, kd, ks, Vn)
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
        %triangle = triangle';
        [sortedV, sortedIndex] = sortrows(triangle); %sort triangle points by y
        ka = ka(sortedIndex, :);                     %rearrabnge colors
        kd = kd(sortedIndex, :); 
        ks = ks(sortedIndex, :); 
        Vn = Vn(sortedIndex, :); 
        ka_grads = zeros(3,3);
        kd_grads = zeros(3,3);
        ks_grads = zeros(3,3);
        Vn_grads = zeros(3,3);
        invSlope = [0 0 0];
        pairs = [1 2; 1 3; 2 3];
        for i = 1:3
            ka_grads(i,:) = (ka(pairs(i,1),:) - ka(pairs(i,2),:))/(sortedV(pairs(i,1),1)-sortedV(pairs(i,2),1));
            kd_grads(i,:) = (kd(pairs(i,1),:) - kd(pairs(i,2),:))/(sortedV(pairs(i,1),1)-sortedV(pairs(i,2),1));
            ks_grads(i,:) = (ks(pairs(i,1),:) - ks(pairs(i,2),:))/(sortedV(pairs(i,1),1)-sortedV(pairs(i,2),1));
            Vn_grads(i,:) = (Vn(pairs(i,1),:) - Vn(pairs(i,2),:))/(sortedV(pairs(i,1),1)-sortedV(pairs(i,2),1));
            
            invSlope(i) = (sortedV(pairs(i,1),1)-sortedV(pairs(i,2),1))/(sortedV(pairs(i,1),2)-sortedV(pairs(i,2),2));
            invSlope(i) = 1/invSlope(i);
        end
        tr = sortedV;
        slopes = invSlope;
    end
    
    function c = colorInterp(Pc, S, C, ncoeff, Ia, I0, ka1, ka2, kd1, kd2, ks1, ks2, vn1, vn2, a, b, x)
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
            ka = ka1;
            kd = kd1;
            ks = ks1;
            vn = vn1;
            c = zeros(3,1);
            c = c + ambientLight(ka', Ia);
            c = c + diffuseLight(Pc, vn', kd', S', I0);
            c = c + specularLight(Pc, vn', C, ks', ncoeff, S', I0);
            

            %diffuseLight(Pc, vn', kd', S', I0)
            c = c';
            return
        end
        %c = A - ((A-B)./(a-b)).*(a - x(:));
        ka = ka1 - ((ka1-ka2)./(a-b)).*(a - x(:));
        kd = kd1 - ((kd1-kd2)./(a-b)).*(a - x(:));
        ks = ks1 - ((ks1-ks2)./(a-b)).*(a - x(:));
        vn = vn1 - ((vn1-vn2)./(a-b)).*(a - x(:));

        c = zeros(3,length(x));
        c = c + ambientLight(ka', Ia);
        c = c + diffuseLight(Pc, vn', kd', S', I0);
        c = c + specularLight(Pc, vn', C, ks', ncoeff, S', I0);
        %size(c)
        c = c';
    end