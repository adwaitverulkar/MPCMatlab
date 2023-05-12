function dydt = non_linear_dynamics(~,x,u,m,l_r,l_f)

    dydt = [x(4)*cos(x(3)) - x(5)*sin(x(3))
       x(4)*sin(x(3)) + x(5)*cos(x(3))
       x(6)
       u(1)/m
       (x(4)*u(2) + x(7)*u(1)/m)*(l_r/(l_r+l_f))
       (x(4)*u(2) + x(7)*u(1)/m)*(1/(l_r+l_f))
       u(2)];

end