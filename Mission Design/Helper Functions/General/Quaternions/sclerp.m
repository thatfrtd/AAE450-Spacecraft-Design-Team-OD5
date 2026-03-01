function [C] = sclerp(A_UDQ, B_UDQ, N)
%SCLERP Screw linear interpolation
%   Interpolate between two dual quaternions

Transf = qualQuatProduct(dualQuatConjucate(A_UDQ), B_UDQ);
[l, m, theta, d] = screwParameters(Transf);

tau = linspace(0, 1, N);
C = zeros(8, length(tau));

for i = 1:length(tau)
    P = [sin(tau(i)*theta/2)*l cos(tau(i)*theta/2)];
    Q = [tau(i)*d/2*cos(tau(i)*theta/2)*l + sin(tau(i)*theta/2)*m -tau(i)*d/2*sin(tau(i)*theta/2)];
    C(:, i) = qualQuatProduct(A_UDQ, [P Q]);
end

end

function R = quat2rot(Q)
q_0 = Q(4);
q_1 = Q(1);
q_2 = Q(2);
q_3 = Q(3);
R(1,1) = q_0*q_0 + q_1*q_1 - q_2*q_2 - q_3*q_3;
R(1,2) = 2*(q_1*q_2 - q_0*q_3);
R(1,3) = 2*(q_0*q_2 + q_1*q_3);
R(2,1) = 2*(q_0*q_3 + q_1*q_2);
R(2,2) = q_0*q_0 - q_1*q_1 + q_2*q_2 - q_3*q_3;
R(2,3) = 2*(q_2*q_3 - q_0*q_1);
R(3,1) = 2*(q_1*q_3 - q_0*q_2);
R(3,2) = 2*(q_0*q_1 + q_2*q_3);
R(3,3) = q_0*q_0 - q_1*q_1 - q_2*q_2 + q_3*q_3;
end

function PQ = quatProduct(P, Q)
p_0 = P(4);
q_0 = Q(4);
p_r = P(1:3);
q_r = Q(1:3);
scalarPart = p_0*q_0 - dot(p_r,q_r);
vectorPart = p_0*q_r + q_0*p_r + cross(p_r,q_r);
PQ = [vectorPart scalarPart];
end

function quatConjucate = quatConjucate(Q)
q_0 = Q(4);
q_r = Q(1:3);
quatConjucate = [-q_r q_0];
end

function dualQuatConjucate = dualQuatConjucate(A)
P = A(1:4);
Q = A(5:8);
dualQuatConjucate = [quatConjucate(P) quatConjucate(Q)];
end

function AB = qualQuatProduct(A, B)
P_A = A(1:4);
Q_A = A(5:8);
P_B = B(1:4);
Q_B = B(5:8);
AB = [quatProduct(P_A,P_B) quatProduct(P_A,Q_B)+quatProduct(Q_A,P_B)];
end

function [u, theta] = quat2AxisAngle(Q)
q_0 = Q(4);
q_r = Q(1:3);
if (norm(q_r) <= 1e-12)
    % Null or full rotation, the angle is 0 (modulo 2*pi) --> singularity: The unit vector u is indeterminate.
    % By convention, set u to the default value [0, 0, 1].
    u = [0, 0, 1];
    theta = 0;
else
    u = q_r/norm(q_r);
    theta = 2*atan2(norm(q_r), q_0);
end
end

function [R, p] = dualQuat2transformation(A)
P = A(1:4);
Q = A(5:8);
R = quat2rot(P);
p_Q = 2 * quatProduct(Q, quatConjucate(P));
p = p_Q(1:3);
end

function [l, m, theta, d] = screwParameters(A)
P = A(1:4);
Q = A(5:8);
[l, theta] = quat2AxisAngle(P);
if theta == 0 || theta == pi
    disp("Screw axis is at infinity!")
else
    [~, t] = dualQuat2transformation(A);
    d = dot(t,l);
    m = 1/2 * (cross(t,l) + (t - d * l) * cot(theta/2));
end
end
