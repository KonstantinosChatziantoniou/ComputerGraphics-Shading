% Konstantinos Chatziantoniou 8941 konstantic@ece.auth.gr
% Aristotle University of Thessaloniki
% Computer Graphics
% 3rd Assignment - 2020/06/19
function Normals = findVertNormals(R, F)
    % params
    %   R   3xr triangle edges.
    %   F   3xT pairs of edges forming a triangle.
    % returns 
    %   Normals 3xr
    [~,r] = size(R);
    [~,T] = size(F);
    Normals = zeros(3,r);
    for i = 1:T
        tr = R(:,F(:,i));
        vecs = [tr(:,2) - tr(:,1), tr(:, 3) - tr(:,1)];
        Normals(:, F(:,i)) = Normals(:, F(:,i)) + cross(vecs(:,1), vecs(:,2));
        %Normals(:, F(:,i)) = zeros(3,3) + cross(vecs(:,1), vecs(:,2));
        
    end
    Normals = normc(Normals);
end