function dX = my_fun(t,x1dot, x2dot)
dX = zeros(2,1);
u  = x1dot;
w  = x2dot;
A  = 1;
B  = 1;
dX = [w*u^2 - B*u;...
      A - w - w*u^2];
end