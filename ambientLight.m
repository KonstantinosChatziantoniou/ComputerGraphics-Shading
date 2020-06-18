% Konstantinos Chatziantoniou 8941 konstantic@ece.auth.gr
% Aristotle University of Thessaloniki
% Computer Graphics
% 3rd Assignment - 2020/06/19
function I = ambientLight(ka, Ia)
% params
%   ka 3xr  ambient light coefs
%   Ia 3x1  ambient light intensity

    I = ka.*Ia;
    
end