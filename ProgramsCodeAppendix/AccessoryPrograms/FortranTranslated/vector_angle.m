function angle = vector_angle(X,Y)

%Returns angle between vectors in rad
%X,Y vectors 1x3, (can take in Nx3 arrays of vectors, then return angle in Nx1).


%finding vector angle:
n = cross(X,Y);
nn = dot(n,n,2);
if nn == 0
    angle = 0;
else
    n = n./sqrt(nn);
    alpha = cross(n,X);
    beta = X - dot(n,X,2).*n;
    gamma = Y - dot(n,Y,2).*n;
    angle = atan2(dot(alpha,gamma,2),dot(beta,gamma,2));
end

end