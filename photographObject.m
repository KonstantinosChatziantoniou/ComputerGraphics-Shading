% Konstantinos Chatziantoniou 8941 konstantic@ece.auth.gr
% Aristotle University of Thessaloniki
% Computer Graphics
% 3rd Assignment - 2020/06/19
function Im = photographObject(shader, f, C, K, u, bC, M, N, H, W, R, F, S, ka, kd, ks, ncoeff, Ia, I0)
    
    NormalVectors = findVertNormals(R, F);
    [P,D] = projectCameraKu(f, C, K, u, R);
    V = rasterize(P, M, N, H, W);
    
    % Calculate the edge of each triangle
    I = ones(M, N, 3);
    for i = 1:3
        I(:,:,i) = I(:,:,i)*bC(i) ;
    end
    t_depth = sum((D(F)),1)/3;
    %    t_depth = min((D(F)));
    % Sort by tringles by depth
    [ ~, depth_index] = sort(-t_depth);

    if shader == 1
        for j = length(depth_index):-1:1
            i = depth_index(j);
            edge_pairs = F(:, i);
            Vn = NormalVectors(:, edge_pairs);
            triangle = V(edge_pairs,:);
            I = shadeGouraud(triangle, Vn,  sum(R(:,edge_pairs),2)/3, C, S, ka(:,edge_pairs), kd(:,edge_pairs), ks(:,edge_pairs), ncoeff, Ia, I0, I);
        end
    end
    if shader == 2
        for j = length(depth_index):-1:1
            i = depth_index(j);
            edge_pairs = F(:, i);
            Vn = NormalVectors(:, edge_pairs);
            triangle = V(edge_pairs,:);
            I = shadePhong(triangle, Vn,  sum(R(:,edge_pairs),2)/3, C, S, ka(:,edge_pairs), kd(:,edge_pairs), ks(:,edge_pairs), ncoeff, Ia, I0, I);
        end
    end
    Im = I;
end



% Konstantinos Chatziantoniou 8941 konstantic@ece.auth.gr
% Aristotle University of Thessaloniki
% Computer Graphics
% 2nd Assignment - 2020/05/15
function [P,D] = projectCamera(w, cv, cx, cy, p)
% Params:
%   w       scalar, distance of the camera from the lense.
%   cv      vector, position of the camera in the world coordinate system.
%   cx      vector, x axis of the camera
%   cy      vector, y axis of the camera
%   p       matrix 3xN, a set of points
% Returns:
%   P       matrix 2xN, the 2d projection of the points
%   D       vector, depth: the distance of the points from the lense.
% Summary:
%   Transforms the points to the coordinate system of the camera. Then project them 
%   to the x-y plane of the camera.
    
    cz = cross(cx,cy);
    p = systemTransform(p, cx, cy, cz, cv);
    
    x = p(1,:);
    y = p(2,:);
    z = p(3,:);
    xp = (w./z).*x;
    xp = xp(:);
    yp = (w./z).*y;
    yp = yp(:);
    P = [xp, -yp];
    D = z(:);
end

% Konstantinos Chatziantoniou 8941 konstantic@ece.auth.gr
% Aristotle University of Thessaloniki
% Computer Graphics
% 2nd Assignment - 2020/05/15
function [P,D] = projectCameraKu(w, cv, ck ,cu, p)
% Params:
%   w       scalar, distance of the camera from the lense.
%   cv      vector, position of the camera in the world coordinate system.
%   ck      vector, position of the target in the world coordinate system.
%   cu      vector, up vector for the cameras orientation.
%   p       matrix 3xN, a set of points
% Returns:
%   P       matrix 2xN, the 2d projection of the points
%   D       vector, depth: the distance of the points from the lense.
% Summary:
%   Calculates the x,y,z axis of the camera given its position, target and up vector.
%   Then calls projectCamera for the projection.
    cu = -cu;
    ck = cv-ck;
    zc = ck/norm(ck);
    yc = cu - dot(cu, zc)*zc;
    yc = yc/norm(yc);
    xc = cross(yc,zc);
    
    [P,D] = projectCamera(w, cv, xc, yc, p);
end

% Konstantinos Chatziantoniou 8941 konstantic@ece.auth.gr
% Aristotle University of Thessaloniki
% Computer Graphics
% 2nd Assignment - 2020/05/15
function Prast = rasterize(P,M,N,H,W)
% Params:
%   P       matrix 2xN, set of projected points.
%   M       scalar, resolution of the camera for x axis
%   N       scalar, resolution of the camera for y axis
%   H       scalar, size of the camera's lense (x axis)
%   W       scalar, size of the camera's lense (y axis)
% Returns:
%   Pras    matrix 2xN, the set of points assigned to a grid.
Prast = [M/H,N/W].*P + [M/2,N/2];
Prast = floor(Prast);


end


% Konstantinos Chatziantoniou 8941 konstantic@ece.auth.gr
% Aristotle University of Thessaloniki
% Computer Graphics
% 2nd Assignment - 2020/05/15
function cq = affineTransform(cp, R, ct)
% Params:
%   cp      matrix 3xN, a set of points in 3d space
%   R       matrix 3x3, rotation matrix
%   ct      vector, for translation
% Returns
%   cq      matrix 3xN, a set of points after rotation and translation
% Summary:
%   Applies rotation and then translation to a set of points and 
%   returns the result.
cq = R*cp;
cq = cq + ct(:);
end

% Konstantinos Chatziantoniou 8941 konstantic@ece.auth.gr
% Aristotle University of Thessaloniki
% Computer Graphics
% 2nd Assignment - 2020/05/15
function dp = systemTransform(cp, b1, b2, b3, c0)
% Params:
%   cp      matrix 3xN, a set of points
%   b1      vector, the x axis coordinates of the new coordinate system
%   b2      vector, the y axis coordinates of the new coordinate system
%   b3      vector, the z axis coordinates of the new coordinate system
%   c0      the translation of the system from the original
% Returns:
%   dp      matrix 3xN, the coordinates of the points to the new system
% Summary: 
%   Apply the inverse rotation and translation to the points to find their
%   coordinates for the new system. We assume the original system is x = [1 0 0]
%   y=[0 1 0] and z = [0 0 1]
L = [b1(:),b2(:),b3(:)];
L = L';
c0 = -L*c0(:);
dp = affineTransform(cp,L, c0);

end