% Konstantinos Chatziantoniou 8941 konstantic@ece.auth.gr
% Aristotle University of Thessaloniki
% Computer Graphics
% 3rd Assignment - 2020/06/19
function I = diffuseLight(P, N, kd, S, I0)
% params
%   P 3x1   3d points (center of mass of triangle)
%   N 3xr   normal vector of each 3d point
%   kd 3xr  diffuse coefs
%   S 3xn   coords of light sources
%   I0 3xn  intensity of each light source
   
    [~, n] = size(S);
    [~, r] = size(N);
    I = zeros(3,r);
    for i = 1:n
        s = S(:,i);
        sI = I0(:,i);
        %L = repmat(P, 1,r) - s;
        L = s- repmat(P, 1,r);
        cosa = dot(N, normc(L));
        cosa = cosa.*(cosa >= 0);
        facc = 1; %1./sum((L).^2);
        I = I + kd.*sI.*cosa.*facc;
    end
    
end
