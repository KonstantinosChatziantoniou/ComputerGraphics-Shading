% Konstantinos Chatziantoniou 8941 konstantic@ece.auth.gr
% Aristotle University of Thessaloniki
% Computer Graphics
% 3rd Assignment - 2020/06/19
function I = specularLight(P, N, C, ks, ncoeff, S, I0)
    % params
    %   P 3xr   3d points
    %   N 3xr   normal vector of each 3d point
    %   C 3x1   camera coords
    %   ks 3xr  specular coefs
    %   ncoeff  scattering level
    %   S 3xn   coords of light sources
    %   I0 3xn  intensity of each light source

    [~, n] = size(S);
    [~, r] = size(N);
    I = zeros(3,r);
    for i = 1:n
        s = S(:,i);
        sI = I0(:,i);
        L = s - repmat(P, 1,r);
        %L = repmat(P, 1,r) - s;
        %V = repmat(P, 1,r) - C;
        V = C - repmat(P, 1,r);
        R = 2*N.*dot(N,normc(L)) - normc(L);
        cosab = dot(normc(R), normc(V));
        cosab = cosab.*(cosab>0);
        cosab = cosab.^ncoeff;
        facc = 1; %1./sum((L).^2)
        I = I + ks.*sI.*cosab.*facc;
    end
end