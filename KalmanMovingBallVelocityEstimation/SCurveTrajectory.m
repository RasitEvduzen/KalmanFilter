function [p,pd,pdd] = SCurveTrajectory(qi, qf, t, V,TS)

tf = max(t);
t = (0:TS:t)';
V = abs(V) * sign(qf-qi);
if abs(V) < abs(qf-qi)/tf
    error('V too small');
elseif abs(V) > 2*abs(qf-qi)/tf
    error('V too big');
end


if qi == qf
    p = ones(size(t)) * qi;
    pd = zeros(size(t));
    pdd = zeros(size(t));
    return
end

tb = (qi - qf + V*tf)/V;
a = V/tb;
p = zeros(length(t), 1);
pd = p;
pdd = p;

for i = 1:length(t)
    
    if t(i) <= tb
        % initial blend
        p(i) = qi + a/2*t(i)^2;
        pd(i) = a*t(i);
        pdd(i) = a;
    elseif t(i) <= (tf-tb)
        % linear motion
        p(i) = (qf+qi-V*tf)/2 + V*t(i);
        pd(i) = V;
        pdd(i) = 0;
    else
        % final blend
        p(i) = qf - a/2*tf^2 + a*tf*t(i) - a/2*t(i)^2;
        pd(i) = a*tf - a*t(i);
        pdd(i) = -a;
    end
end

end