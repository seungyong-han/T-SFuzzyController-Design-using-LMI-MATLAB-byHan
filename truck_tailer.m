function dx=truck_tailer(x,u)
global v tbar L l t0
dx(1,1) =  -((v*tbar)/(L*t0))*x(1) + ((v*tbar)/(l*t0))*u;
dx(2,1) = ((v*tbar)/(L*t0))*x(1);
dx(3,1) = ((v*tbar)/t0)*sin(x(2)+((v*tbar)/(2*L*t0))*x(1));


