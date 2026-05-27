%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   A solution branch code for the standard non-dimensional blood flow model Riemann problem.
%  
%    Bethel fecit, AD 2025.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, clear the workspace and close any open figures
clearvars
clc
close all


% Set precision
format short

% Figure
figure 
box on

% Set initial conditions for the Riemann problem
u1 = 0.2;
A1 = 0.8;
u2 = 0.4;
A2bar= 1;
k2=1.1;
m=10;
n=0;
dur2=0.005;
A2range= 0.05:dur2:2;


%initial guess for the resonant solutions
A11_guess = 1.05;
A12_guess = A11_guess;
guess_Ro0oS = [0.0695    0.4126    0.7675    0.7935];
guess_Ro0oR = [0.0905    0.0993    0.6879    0.7094];
guess_RoSS = [0.4126    0.7935];
guess_RoSR = [0.4126    0.7935];

% Define the loop for A2
for  An = 1:length(A2range)
        A2 = A2range(An);
        disp('A2:');
        disp(A2);

 c2=sqrt(k2*(m*(A2/A2bar).^m+n*(A2/A2bar).^(-n)));

  if u2<c2 && u2>-c2  % Consider only subcritical regime

    [x, exitflag] = classical_sing_gen(u1,u2,k2,m,A2bar,A1,A2,n);
    [exitflag_RoSS,guess_RoSS] = resonance_RoSS_sing_gen(u1,u2,k2,m,n,A2bar,A1,A2,guess_RoSS,A12_guess);
    [exitflag_RoSR,guess_RoSR] = resonance_RoSR_sing_gen(u1,u2,k2,m,n,A2bar,A1,A2,guess_RoSR,A12_guess);
    [exitflag_Ro0oS,guess_Ro0oS] = resonance_Ro0oS_sing_gen(u1,u2,k2,m,n,A2bar,A1,A2,guess_Ro0oS,A11_guess);
    [exitflag_Ro0oR,guess_Ro0oR] = resonance_Ro0oR_sing_gen(u1,u2,k2,m,n,A2bar,A1,A2,guess_Ro0oR,A11_guess);
           
      %Plot A_11
 figure(1)
    box on
    subplot 211
      
        if  exitflag==1 && x(3)>0 && x(1)>0 
                
            if x(1) <= A1 && x(3) <= A2  % RoR 
                hold on
                plot(A2,x(1), 'b.');    % Plots blue dots on parameter space graph
        
            elseif x(1) <= A1 && x(3) > A2 % RoS
                hold on
                plot(A2,x(1), 'g.');    % Plots green dots on parameter space graph
        
            elseif x(1) > A1 && x(3) <= A2  % SoR
                hold on
                plot(A2,x(1), 'c.');   % Plots cyan dots on parameter space graph
        
            elseif x(1) > A1 && x(3) > A2 % SoS
                hold on
                plot(A2,x(1), 'r.');   % Plots red dots on parameter space graph
           
            end
       elseif exitflag_RoSS == 1 && guess_RoSS(1)>0 && guess_RoSS(2)>0 
            hold on
            plot(A2,guess_RoSS(3), '.',"MarkerEdgeColor", [0.8500 0.3250 0.0980], "MarkerFaceColor",[0.8500 0.3250 0.0980]); % Plots orange dots on parameter space graph
        
        elseif exitflag_RoSR == 1 && guess_RoSR(1)>0 && guess_RoSR(2)>0 
            hold on
            plot(A2,guess_RoSR(3), 'm.'); % Plot magenta dots on bifurcation graph

        elseif exitflag_Ro0oS == 1 && guess_Ro0oS(1)>0 && guess_Ro0oS(2)>0 && guess_Ro0oS(3)>0 && guess_Ro0oS(4)>0
            hold on
            plot(A2,guess_Ro0oS(5), '.',"MarkerEdgeColor", "#7e1e9c", "MarkerFaceColor","#7e1e9c"); % Plots purple dots on parameter space graph
       
        elseif exitflag_Ro0oR == 1 && guess_Ro0oR(1)>0 && guess_Ro0oR(2)>0 && guess_Ro0oR(3)>0 && guess_Ro0oR(4)>0
            hold on
            plot(A2,guess_Ro0oR(5), 'y.'); % Plots yellow dots on parameter space graph

        else                    % No solution
            hold on                     
            plot(A2,u2, 'k.');   % Plots black cross on parameter space graph
       end
      
    xlabel('$A_2$','interpreter','latex','FontSize',11)
    ylabel('$A_{11}$','interpreter','latex','FontSize',11)


 ylim([0 2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot2 A_22

    figure(1)
    box on
subplot 212
        if  exitflag==1 && x(3)>0 && x(1)>0 
                
            if x(1) <= A1 && x(3) <= A2  % RoR 
                hold on
                plot(A2,x(3), 'b.');    % Plots blue dots on parameter space graph
        
            elseif x(1) <= A1 && x(3) > A2 % RoS
                hold on
                plot(A2,x(3), 'g.');    % Plots green dots on parameter space graph
        
            elseif x(1) > A1 && x(3) <= A2  % SoR
                hold on
                plot(A2,x(3), 'c.');   % Plots cyan dots on parameter space graph
        
            elseif x(1) > A1 && x(3) > A2 % SoS
                hold on
                plot(A2,x(3), 'r.');   % Plots red dots on parameter space graph
           
            end
       elseif exitflag_RoSS == 1 && guess_RoSS(1)>0 && guess_RoSS(2)>0 
            hold on
            plot(A2,guess_RoSS(2), '.',"MarkerEdgeColor", [0.8500 0.3250 0.0980], "MarkerFaceColor",[0.8500 0.3250 0.0980]); % Plots orange dots on parameter space graph
        
        elseif exitflag_RoSR == 1 && guess_RoSR(1)>0 && guess_RoSR(2)>0 
            hold on
            plot(A2,guess_RoSR(2), 'm.'); % Plot magenta dots on bifurcation graph

        elseif exitflag_Ro0oS == 1 && guess_Ro0oS(1)>0 && guess_Ro0oS(2)>0 && guess_Ro0oS(3)>0 && guess_Ro0oS(4)>0
            hold on
            plot(A2,guess_Ro0oS(4), '.',"MarkerEdgeColor", "#7e1e9c", "MarkerFaceColor","#7e1e9c"); % Plots purple dots on parameter space graph
       
        elseif exitflag_Ro0oR == 1 && guess_Ro0oR(1)>0 && guess_Ro0oR(2)>0 && guess_Ro0oR(3)>0 && guess_Ro0oR(4)>0
            hold on
            plot(A2,guess_Ro0oR(4), 'y.'); % Plots yellow dots on parameter space graph

        else                    % No solution
            hold on                     
            plot(A2,u2, 'k.');   % Plots black cross on parameter space graph
       end 
 
    xlabel('$A_2$','interpreter','latex','FontSize',11)
    ylabel('$A_{22}$','interpreter','latex','FontSize',11)

ylim([0 2])

   end
 end %A2
 
  