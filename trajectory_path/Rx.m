function T=Rx(qx)
    T = [1   0        0     0;
         0 cos(qx) -sin(qx) 0;
         0 sin(qx)  cos(qx) 0;
         0   0        0     1];
end