% q=sym('q', [6,1]);
% assume(q, 'real')
% q;
% B = [-eye(3), eye(3)];
% % B =sym('B', [3, 6]);
% syms k l0
% V = ((sqrt((q.')*(B.')*B*q) - l0)^2)%*k/2;
% 
% B
% 
% sqrt((q.')*(B.')*B*q)
% 
% ccode(simplify(gradient(V, q)))
% ccode(gradient(V, q))
% 
% ccode(simplify(hessian(V, q)))
% ccode(hessian(V, q))
% 
% ccode(gradient(gradient(V, q), q))

%%%%%%%%%%%%%%%%%%%%%%%
q0 = sym('q0', [3,1]);
q1 = sym('q1', [3,1]);
assume(q0, 'real')
assume(q1, 'real')
syms k l0
 V = ((norm(q0-q1) - l0)^2)*k/2;

ccode((simplify(gradient(V, q0))))
ccode((simplify(gradient(V, q1))))

ccode((simplify(hessian(V, q0))))
ccode(hessian(V, q0))