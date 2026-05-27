%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   A parameter space code for the standard non-dimensional blood flow model Riemann problem using general tube law.
%
%   Bethel fecit, AD 2025.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% First, clear the workspace and close any open figures
clc
close all
clear

% Set precision
format short

% Figure
figure(1)
box on

% Set initial conditions for the Riemann problem
A1 = 0.8;
u1 = 0.2;
A2bar= 1;
k2=1.1;
m=10; % Go to classical_sinf_gen.m and choose the appropriate initial guess for the classical case
n=1.5;

% Initial guesses for resonant cases
A11_guess = 1.05;
A12_guess = A11_guess;
guess_Ro0oS = [0.10905    0.0759    0.1674    0.1722];
guess_Ro0oR = [0.2587    0.3765    0.4959    0.5630];
guess_RoSS = [0.0508    0.1486];
guess_RoSR = [0.2741    0.4020];


% Find the solutions for different u2 and A2, and plot them 
% First set the ranges of u2 and A2
dur1 = 0.025;
u2range =-2.5:dur1:2.5;

dur2=0.01;
A2range=0.06:dur2:2;

% Define the loops for u2 and A2
for un = 1:length(u2range)
    for An = 1:length(A2range)
        u2 = u2range(un);
        A2 = A2range(An);
        disp('A2:');
        disp(A2);

    c2=sqrt(k2*(m*(A2/A2bar).^m+n*(A2/A2bar).^(-n)));
    
 if u2<c2 && u2>-c2  % Consider only subcritical regime

          
    % Decide what type of solution

  figure(1)
       [x, exitflag] = classical_sing_gen(u1,u2,k2,m,A2bar,A1,A2,n);
     
     if  exitflag==1 && x(3)>0 && x(1)>0 
                
            if x(1) <= A1 && x(3) <= A2  % RoR 
                hold on
                plot(A2,u2, 'b.');    % Plots blue dots on parameter space graph
        
            elseif x(1) <= A1 && x(3) > A2 % RoS
                hold on
                plot(A2,u2, 'g.');    % Plots green dots on parameter space graph
        
            elseif x(1) > A1 && x(3) <= A2  % SoR
                hold on
                plot(A2,u2, 'c.');   % Plots cyan dots on parameter space graph
        
            elseif x(1) > A1 && x(3) > A2 % SoS
                hold on
                plot(A2,u2, 'r.');   % Plots red dots on parameter space graph
           
            end
     else
        [exitflag_RoSS,guess_RoSS] = resonance_RoSS_sing_gen(u1,u2,k2,m,n,A2bar,A1,A2,guess_RoSS,A12_guess);
        [exitflag_RoSR,guess_RoSR] = resonance_RoSR_sing_gen(u1,u2,k2,m,n,A2bar,A1,A2,guess_RoSR,A12_guess);
        [exitflag_Ro0oS,guess_Ro0oS] = resonance_Ro0oS_sing_gen(u1,u2,k2,m,n,A2bar,A1,A2,guess_Ro0oS,A11_guess);
        [exitflag_Ro0oR,guess_Ro0oR] = resonance_Ro0oR_sing_gen(u1,u2,k2,m,n,A2bar,A1,A2,guess_Ro0oR,A11_guess);
    
        if exitflag_RoSS == 1 && guess_RoSS(1)>0 && guess_RoSS(2)>0 
            hold on
            plot(A2,u2, '.',"MarkerEdgeColor", [0.8500 0.3250 0.0980], "MarkerFaceColor",[0.8500 0.3250 0.0980]); % Plots orange dots on parameter space graph
        
        elseif exitflag_RoSR == 1 && guess_RoSR(1)>0 && guess_RoSR(2)>0 
            hold on
            plot(A2,u2, 'm.'); % Plot magenta dots on bifurcation graph

        elseif exitflag_Ro0oS == 1 && guess_Ro0oS(1)>0 && guess_Ro0oS(2)>0 && guess_Ro0oS(3)>0 && guess_Ro0oS(4)>0
            hold on
            plot(A2,u2, '.',"MarkerEdgeColor", "#7e1e9c", "MarkerFaceColor","#7e1e9c"); % Plots purple dots on parameter space graph
       
        elseif exitflag_Ro0oR == 1 && guess_Ro0oR(1)>0 && guess_Ro0oR(2)>0 && guess_Ro0oR(3)>0 && guess_Ro0oR(4)>0
            hold on
            plot(A2,u2, 'y.'); % Plots yellow dots on parameter space graph

        else                    % No solution
            hold on                     
            plot(A2,u2, 'k.');   % Plots black cross on parameter space graph
       end
      end
 end
    xlabel('$A_2$','interpreter','latex','FontSize',11)
    ylabel('$u_2$','interpreter','latex','FontSize',11)
   
     end %A2

end % u2

% Critical boundary
c2=sqrt(k2*(m*(A2range/A2bar).^m+n*(A2range/A2bar).^(-n)));
hold on
plot(A2range,c2,'k')
ylim([-2.5 2.5])
plot(A2range,-c2,'k')
