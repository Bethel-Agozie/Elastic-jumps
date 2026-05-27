function [x, exitflag]= classical_sing_gen(u1,u2,k2,m1,m2,A2bar,A1,A2,n1,n2)

% Define the system of nonlinear equations
fun = @(x) [
    x(2) - u1 + piecewise_l_gen(x,A1,m1,n1);
    x(1)*x(2) - x(3)*x(4);
    0.5*(x(2)^2 - x(4)^2) + x(1)^m1-x(1)^(-n1)-k2*((x(3)/A2bar)^m2-(x(3)/A2bar)^(-n2));
    x(4)- u2 - piecewise_r_gen(x,A2,k2,m2,n2,A2bar)
];

% Initial guess

x0 = [1; 1; 0.8; 1.2]; 

% Solve the system of equations using fsolve
options = optimoptions('fsolve', 'Display', 'iter', 'MaxIterations', 10000, 'FunctionTolerance', 1e-8);
[x,fval,exitflag] = fsolve(fun, x0, options);

% Display the result
disp('Solution:');
disp(x);
disp(exitflag);

if exitflag==1 && x(1)>0 && x(3)>0 

    if x(1) <= A1 && x(3) <= A2
        disp('RoR')
        
    elseif x(1) <= A1 && x(3) > A2  
           disp('RoS')
                       
    elseif x(1) > A1 && x(3) <= A2
        disp('SoR')
                
    elseif x(1) > A1 && x(3) > A2
        disp('SoS')
              
    end
else
    disp('Not solved')
end


