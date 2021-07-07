function R = vector_rotate(X, A, Phi) %X input vector, A rotation axis, Phi angle.
%vecotr input format Nx3

norm = A./sqrt(dot(A,A,2));

%R is the rotated vector
R = norm.*(dot(norm,X,2)) + cross(norm, X).*sin(Phi) + (X - norm.*dot(norm,X,2)).*cos(Phi);

end