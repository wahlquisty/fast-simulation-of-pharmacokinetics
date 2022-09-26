function y = simulate_lsim(PK,uinf,ubol,t,youts)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
nt = length(t);
y = zeros(nt,1);
x0 = [0,0,0]';
x = x0;

for i = 1:nt-1
    x = x + PK.B * ubol(i); % add bolus, it this correct?!
    [Y,T,X] = lsim(PK,[uinf(i) uinf(i+1)],[t(i) t(i+1)],x);
    y(i) = Y(1);
    x = X(end,:)';
end
% last datapoint
x = x + PK.B * ubol(nt); % add bolus, it this correct?!
[Y,T,X] = lsim(PK,[uinf(nt-1) uinf(nt)],[t(nt-1) t(nt)],x);
y(nt) = Y(1);
x = X(end,:)';
    
y = y(youts);
end

