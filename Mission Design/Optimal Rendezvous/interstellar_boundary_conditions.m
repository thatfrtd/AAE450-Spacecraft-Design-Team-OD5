function [x_f, x_array] = interstellar_boundary_conditions(q_0, w_0, t_k, MOI, d_vec) 

q_dvec_body = point_at_vec(-d_vec);

odeop = odeset(RelTol=1e-13,AbsTol=1e-13);

[~,x_array] = ode45(@(t,x)ACS(x,MOI),t_k,[q_0;w_0], odeop);
x = x_array(end,:);
x_array = x_array';

q_f = x(end,1:4)';

w_f = x(end,5:7)';

q_dvec_body_f = q_mul(q_f,q_dvec_body);

w_quac = quat_rot(q_conj(q_dvec_body_f),w_f);

[d_fin] = quat_rot(q_f,d_vec);

Iner_Ve = quat_rot(q_f, cross(w_f,d_vec)); 

x_f = [d_fin;Iner_Ve;q_dvec_body_f;w_quac];

end


%%Helper Func
function x_dot = ACS(x,I) 
% Derrivative of state
q = x(1:4);
w = x(5:7);
tau = [0;0;0];
q_dot = QKE(q,w);
w_dot = Euler_BB(w,I,tau);

x_dot = [q_dot;w_dot];

end

function [quat] = point_at_vec(T_vec)
    T_vec = T_vec ./ vecnorm(T_vec);
    w = 1 + T_vec(1, :);
    xyz = cross([ones([1, size(T_vec, 2)]); zeros([2, size(T_vec, 2)])], T_vec); % cross b_x with T_vec
    quat = [xyz; w];
    quat = quat ./ vecnorm(quat);
end

function [w_dot] = Euler_BB(w, I, M)
    w_dot = (M([1; 2; 3]) + (I([2; 3; 1]) - I([3; 1; 2])) .* w([2; 3; 1]) .* w([3; 1; 2])) ./ I([1; 2; 3]);
end

function [qdot] = QKE(q, w)
    qdot = 1 / 2 * q_mul(q, [w; 0]);
end