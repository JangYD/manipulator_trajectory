%robot dh parameter
function dh_parameter = dh_parameter(a, alpha, d, theta)
    theta = theta(:);
    d = d(:);
    alpha = alpha(:);
    a = a(:);
    z = {a,alpha,d,theta};
    dh_parameter = cat(2,z{:});
end