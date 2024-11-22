
%<USER>
% gen.coordinates
syms th1 th2 real %gc
q = [th1;th2];

syms om1 om2 real %gv
u = [om1;om2];

syms al1 al2 real %gvdot
udot = [al1;al2];

% constants
syms m1 m2 real
syms l1 l2 real
syms g real

% define Lagrangian
r1 = [l1*sin(th1);-l1*cos(th1)];
r2 = r1 + [l2*sin(th2);-l2*cos(th2)];

v1 = [l1*cos(th1)*om1;l1*sin(th1)*om1];
v2 = v1 + [l2*cos(th2)*om2;l2*sin(th2)*om2];

L  = 1/2* ( m1*transpose(v1)*v1 + m2*transpose(v2)*v2) - g*(m1*r1(2) + m2*r2(2));
%<\USER>

el_pot       = transpose(simplify(jacobian(L,q)));
el_kin_no_dt = transpose(simplify(jacobian(L,u)));

% <USER> need t-representation for time derivatives
syms om1t(t) om2t(t) th1t(t) th2t(t);
% function subs. for time derivative
forw_orig = {om1 ,om2 ,th1 ,th2 };
forw_time = {om1t,om2t,th1t,th2t};
% for reverse subs. after time der.
rev_time = {diff(om1t(t),t),diff(om2t(t),t),om1t(t),diff(th1t(t),t),om2t(t),diff(th2t(t),t),th1t(t),th2t(t)};
rev_orig = {al1            ,            al2,    om1,            om1,    om2,            om2,    th1,    th2};
% </USER>

el_kin_no_dt_t = subs(el_kin_no_dt,forw_orig,forw_time);
el_kin_t = simplify(diff(el_kin_no_dt_t,t));

el_kin = subs(el_kin_t,rev_time,rev_orig);
el_kin = simplify(el_kin);

el_full = simplify(el_kin - el_pot);

el_full_mass = simplify(transpose(jacobian(el_full,udot)));
el_full_bias = simplify(el_full - el_full_mass*udot);

disp("Mass Matrix:");
disp(el_full_mass);
disp("Bias Force");
disp(el_full_bias);