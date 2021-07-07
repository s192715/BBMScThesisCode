function impact = impact_parameter(r_leo,r_gps)

%takes in position vectors in ECEF Nx3

%length of LEO and GPS vectors
r0 = sqrt(dot(r_leo,r_leo,2)); %Format: Nx1
r1 = sqrt(dot(r_gps,r_gps,2));

omega = pi - vector_angle(r_leo,r_gps); %in rad

talpha = r1.*sin(omega)./(r0 + r1.*cos(omega));

impact = r0.*talpha./sqrt(1 + talpha.^2);

end