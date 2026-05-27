%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   A wave profile code for the standard non-dimensional blood flow model Riemann using general tube law.
%
%   Bethel fecit, AD 2025.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, clear the workspace and close any open figures
clearvars
clc
close all

% Set precision
format short

% Set initial conditions for the Riemann problem
u1 = 0.2;
A1 = 0.8;

u2 = 0.4;
A2 = 0.88;

A2bar= 1;
k2=1.1;
m=10;
n=0;
t=0.8;

% Initial guess for the resonant solutions
A11_guess = 1.05;
A12_guess = A11_guess;
guess_Ro0oS = [0.0695    0.4126    0.7675    0.7935];
guess_Ro0oR = [0.0905    0.0993    0.6879    0.7094];
guess_RoSS = [0.4126    0.7935];
guess_RoSR = [0.4126    0.7935];


 c1=sqrt(m*A1.^m+n*A1.^(-n));
 c2=sqrt(k2*(m*(A2/A2bar).^m+n*(A2/A2bar).^(-n)));
 
 if u2<c2 && u2>-c2  % Consider only subcritical regime 

 % solve the equations and decide what type of solution
    [exitflag_RoSS,guess_RoSS] = resonance_RoSS_sing_gen(u1,u2,k2,m,n,A2bar,A1,A2,guess_RoSS,A12_guess);
    [exitflag_RoSR,guess_RoSR] = resonance_RoSR_sing_gen(u1,u2,k2,m,n,A2bar,A1,A2,guess_RoSR,A12_guess);
    [exitflag_Ro0oS,guess_Ro0oS] = resonance_Ro0oS_sing_gen(u1,u2,k2,m,n,A2bar,A1,A2,guess_Ro0oS,A11_guess);
    [exitflag_Ro0oR,guess_Ro0oR] = resonance_Ro0oR_sing_gen(u1,u2,k2,m,n,A2bar,A1,A2,guess_Ro0oR,A11_guess);
    
   [x, exitflag] = classical_sing_gen(u1,u2,k2,m,A2bar,A1,A2,n);

%  Compute celerities, shock speeds, decide what type of solution and plot
 
% Plot cross-sectional area, A 
    figure('Position', [234, 150, 350, 240]) % Set the pixel size of the figure
    box on
    if exitflag_Ro0oS == 1 && guess_Ro0oS(1)>0 && guess_Ro0oS(2)>0 && guess_Ro0oS(3)>0 && guess_Ro0oS(4)>0
            c11=sqrt((m*guess_Ro0oS(5).^m+n*guess_Ro0oS(5).^(-n)));
            ReS=(guess_Ro0oS(2)*guess_Ro0oS(7)-guess_Ro0oS(3)*guess_Ro0oS(8))/(guess_Ro0oS(2)-guess_Ro0oS(3));
            rS=(A2*u2-guess_Ro0oS(4)*guess_Ro0oS(9))/(A2-guess_Ro0oS(4));

            % distances 
            x1=(u1-c1)*t; 
            x11 = (guess_Ro0oS(6)-c11)*t;
            x12=x11;
            xReS=ReS*t; 
            x21=xReS;
            x22=x21;
            x2 = rS*t;

               % plots
        X=linspace(x1-1,x1,20); a=linspace(A1,A1,20); plot(X, a, 'color', "#f97306",LineWidth=1);
        hold on
        X=linspace(x1,x11,20); a=linspace(A1,guess_Ro0oS(5),20);plot(X, a, 'b',LineWidth=1);
        hold on
        X=linspace(x11,x12,20); a=linspace(guess_Ro0oS(5), guess_Ro0oS(2),20); plot(X, a, 'mo',LineWidth=1.5);
        hold on
        X=linspace(x12,xReS,20); a=linspace(guess_Ro0oS(2), guess_Ro0oS(2),20); plot(X, a, 'ko',LineWidth=1.5);
        hold on
        X=linspace(xReS,xReS,20); a=linspace(guess_Ro0oS(2), guess_Ro0oS(3),20); plot(X, a, 'y*',LineWidth=1.5);
        hold on
        X=linspace(xReS,x21,20); a=linspace(guess_Ro0oS(3), guess_Ro0oS(3),20); plot(X, a, 'go',LineWidth=1.5);
        hold on
        X=linspace(x21,x22,20); a=linspace(guess_Ro0oS(3), guess_Ro0oS(4),20); plot(X, a, ':','color', "#7e1e9c",LineWidth=2);
        hold on
        X=linspace(x22,x2,20); a=linspace(guess_Ro0oS(4), guess_Ro0oS(4),20); plot(X, a, 'c',LineWidth=1);
        hold on
        X=linspace(x2,x2,20); a=linspace(guess_Ro0oS(4), A2,20); plot(X, a, 'r',LineWidth=1);
        hold on
        X=linspace(x2,x2+1,20); a=linspace(A2, A2,20); plot(X, a, 'color', "#6e750e",LineWidth=1);
        legend({'$A_1$','1-$R$', '3-$o_l$','$A_{12}$','$S_{=0}$','$A_{21}$','3-$o_r$',...
                    '$A_{22}$','2-$S$', '$A_2$'},...
                    'interpreter','latex','FontSize',12,'Location','northeastoutside','Orientation','vertical');

     elseif exitflag_Ro0oR == 1 && guess_Ro0oR(1)>0 && guess_Ro0oR(2)>0 && guess_Ro0oR(3)>0 && guess_Ro0oR(4)>0
        c11=sqrt((m*guess_Ro0oR(5).^m+n*guess_Ro0oR(5).^(-n)));
        c22=sqrt(k2*(m*(guess_Ro0oR(4)/A2bar).^m+n*(guess_Ro0oR(4)/A2bar).^(-n)));
        ReS=(guess_Ro0oR(2)*guess_Ro0oR(7)-guess_Ro0oR(3)*guess_Ro0oR(8))/(guess_Ro0oR(2)-guess_Ro0oR(3));
          % distances 
        x1=(u1-c1)*t; 
        x11 = (guess_Ro0oR(6)-c11)*t;
        x12=x11;
        xReS=ReS*t;
        x21=xReS;
        x22 = (guess_Ro0oR(9) + c22)*t;
        x2 = (u2 + c2)*t;

            % plots
    X=linspace(x1-1,x1,20); a=linspace(A1,A1,20); plot(X, a, 'color', "#f97306",LineWidth=1.5);
    hold on
    X=linspace(x1,x11,20); a=linspace(A1,guess_Ro0oR(5),20);plot(X, a, 'b',LineWidth=1.5);
    hold on
    X=linspace(x11,x12,20); a=linspace(guess_Ro0oR(5), guess_Ro0oR(2),20); plot(X, a, 'mo',LineWidth=1.5);
    hold on
    X=linspace(x12,xReS,20); a=linspace(guess_Ro0oR(2), guess_Ro0oR(2),20); plot(X, a, 'ko',LineWidth=1.5);
    hold on
    X=linspace(xReS,xReS,20); a=linspace(guess_Ro0oR(2), guess_Ro0oR(3),20); plot(X, a, 'y*',LineWidth=1.5);
    hold on
    X=linspace(xReS,x21,20); a=linspace(guess_Ro0oR(3), guess_Ro0oR(3),20); plot(X, a, 'go',LineWidth=1.5);
    hold on
    X=linspace(x21,x21,20); a=linspace(guess_Ro0oR(3), guess_Ro0oR(4),20); plot(X, a,':','color', "#7e1e9c",LineWidth=2)
    hold on
    X=linspace(x21,x22,20); a=linspace(guess_Ro0oR(4), guess_Ro0oR(4),20); plot(X, a, 'c',LineWidth=1.5);
    hold on
    X=linspace(x22,x2,20); a=linspace(guess_Ro0oR(4), A2,20); plot(X, a, 'r',LineWidth=1.5);
    hold on
    X=linspace(x2,x2+1,20); a=linspace(A2, A2,20); plot(X, a, 'color', "#6e750e",LineWidth=1.5);

            legend({'$A_1$','1-$R$', '3-$o_l$','$A_{12}$','$S_{=0}$','$A_{21}$','3-$o_r$',...
                '$A_{22}$','2-$R$', '$A_2$'},...
                'interpreter','latex','FontSize',12,'Location','northeastoutside','Orientation','vertical');

    elseif exitflag_RoSS == 1 && guess_RoSS(1)>0 && guess_RoSS(2)>0  
        c12=sqrt((m*guess_RoSS(4).^m+n*guess_RoSS(4).^(-n)));
        ReS=(guess_RoSS(2)*guess_RoSS(3)-guess_RoSS(1)*guess_RoSS(6))/(guess_RoSS(2)-guess_RoSS(1));
        rS=(A2*u2-guess_RoSS(2)*guess_RoSS(3))/(A2-guess_RoSS(2));

       % distances 
        x1=(u1-c1)*t; 
        x12 = (guess_RoSS(5)-c12)*t;
        x21=x12;
        xReS = ReS*t;
        x2 = rS*t;
       
       % plots
        X=linspace(x1-1,x1,20); a=linspace(A1,A1,20); plot(X, a, 'color', "#f97306",LineWidth=1);
        hold on
        X=linspace(x1,x12,20); a=linspace(A1,guess_RoSS(4),20);plot(X, a, 'b',LineWidth=1);
        hold on
        X=linspace(x12,x21,20); a=linspace(guess_RoSS(4),guess_RoSS(1),20);plot(X, a, 'm:',LineWidth=1);
        hold on
        X=linspace(x21,xReS,20); a=linspace(guess_RoSS(1), guess_RoSS(1),20); plot(X, a, 'g',LineWidth=1);
        hold on
        X=linspace(xReS,xReS,20); a=linspace(guess_RoSS(1), guess_RoSS(2),20); plot(X, a, 'color', "#7e1e9c",LineWidth=1);
        hold on
        X=linspace(xReS,x2,20); a=linspace(guess_RoSS(2), guess_RoSS(2),20); plot(X, a, 'c',LineWidth=1);
        hold on
        X=linspace(x2,x2,20); a=linspace(guess_RoSS(2), A2,20); plot(X, a, 'r',LineWidth=1);
        hold on
        X=linspace(x2,x2+1,20); a=linspace(A2, A2,20); plot(X, a, 'color', "#6e750e",LineWidth=1);
            legend({'$A_1$','1-$R$', '3-$o$','$A_{21}$','$S_{\neq 0}$',...
                '$A_{22}$','2-$S$', '$A_2$'},...
                'interpreter','latex','FontSize',12,'Location','northeastoutside','Orientation','vertical');

     elseif exitflag_RoSR == 1 && guess_RoSR(1)>0 && guess_RoSR(2)>0
        c12=sqrt(m*guess_RoSR(4).^m+n*guess_RoSR(4).^(-n));
        c22=sqrt(k2*(m*(guess_RoSR(2)/A2bar).^m+n*(guess_RoSR(2)/A2bar).^(-n)));
        ReS=(guess_RoSR(2)*guess_RoSR(3)-guess_RoSR(1)*guess_RoSR(6))/(guess_RoSR(2)-guess_RoSR(1));
        
        % distances 
        x1=(u1-c1)*t; 
        x12 = (guess_RoSR(5)-c12)*t;
        xReS = ReS*t;
        x22 = (guess_RoSR(3) + c22)*t; 
        x2=(u2+c2)*t;
       
        % plots
            X=linspace(x1-1,x1,20); a=linspace(A1,A1,20); plot(X, a, 'color', "#f97306",LineWidth=1);
            hold on
            X=linspace(x1,x12,20); a=linspace(A1,guess_RoSR(4),20);plot(X, a, 'b',LineWidth=1);
             hold on
            X=linspace(x12,x12,20); a=linspace(guess_RoSR(4),guess_RoSR(1),20);plot(X, a, 'm:',LineWidth=1);
            hold on
            X=linspace(x12,xReS,20); a=linspace(guess_RoSR(1), guess_RoSR(1),20); plot(X, a, 'g',LineWidth=1);
            hold on
            X=linspace(xReS,xReS,20); a=linspace(guess_RoSR(1), guess_RoSR(2),20); plot(X, a,  'color', "#7e1e9c",LineWidth=1);
            hold on
            X=linspace(xReS,x22,20); a=linspace(guess_RoSR(2), guess_RoSR(2),20); plot(X, a, 'c',LineWidth=1);
            hold on
            X=linspace(x22,x2,20); a=linspace(guess_RoSR(2), A2,20); plot(X, a, 'r',LineWidth=1);
            hold on
            X=linspace(x2,x2+1,20); a=linspace(A2, A2,20); plot(X, a, 'color', "#6e750e",LineWidth=1); 

            legend({'$A_1$','1-$R$', '3-$o$','$A_{21}$','$S_{\neq 0}$',...
                '$A_{22}$','2-$R$', '$A_2$'},...
                'interpreter','latex','FontSize',12,'Location','northeastoutside','Orientation','vertical');
            
    elseif exitflag_So0oR == 1 && guess_So0oR(1)>0 && guess_So0oR(2)>0 && guess_So0oR(3)>0 && guess_So0oR(4)>0
        lS=(A1*u1-guess_So0oR(2)*guess_So0oR(7))/(A1-guess_So0oR(2));
        ReS=(guess_So0oR(3)*guess_So0oR(8)-guess_So0oR(4)*guess_So0oR(9))/(guess_So0oR(3)-guess_So0oR(4));
        c22=sqrt(k2*(m*guess_So0oR(5).^m+n*guess_So0oR(5).^(-n)));
        % distances 
        x1=lS*t; 
        xReS=ReS*t;
        x12=xReS;
        %x11=x12;
        
        x22=(guess_So0oR(6)+c22)*t;
        x21=x22;
        x2=(u2+c2)*t;
                
          % plots
    
        X=linspace(x1-1,x1,20); a=linspace(A1,A1,20); plot(X, a, 'color', "#f97306",LineWidth=1);
        hold on
        X=linspace(x1,x1,20); a=linspace(A1,guess_So0oR(2),20);plot(X, a, 'b',LineWidth=1);
        hold on
        X=linspace(x1,x12,20); a=linspace(guess_So0oR(2),guess_So0oR(2),20);plot(X, a, 'g',LineWidth=1);
        hold on
        X=linspace(x12,x12,20); a=linspace(guess_So0oR(2), guess_So0oR(3),20); plot(X, a, ':','color', "#7e1e9c",LineWidth=2);
        hold on
        X=linspace(x12,xReS,20); a=linspace(guess_So0oR(3), guess_So0oR(3),20); plot(X, a, 'go',LineWidth=1.5);
        hold on
        X=linspace(xReS,xReS,20); a=linspace(guess_So0oR(3), guess_So0oR(4),20); plot(X, a, 'y*',LineWidth=1.5);
        hold on
        X=linspace(xReS,x21,20); a=linspace(guess_So0oR(4), guess_So0oR(4),20); plot(X, a, 'ko',LineWidth=1.5);
        hold on
        X=linspace(x21,x22,20); a=linspace(guess_So0oR(4), A22,20); plot(X, a, 'm',LineWidth=1.5);
        hold on
        X=linspace(x22,x2,20); a=linspace(A22, A2,20); plot(X, a, 'r',LineWidth=1);
        hold on
        X=linspace(x2,x2+1,20); a=linspace(A2, A2,20); plot(X, a, 'color', "#6e750e",LineWidth=1)
        legend({'$A_1$','1-$S$','$A_{11}$','3-$o_l$','$A_{12}$','$S_{=0}$','$A_{21}$','3-$o_r$',...
                    '2-$R$', '$A_2$'},...
                    'interpreter','latex','FontSize',12,'Location','northeastoutside','Orientation','vertical');
    
  elseif exitflag_SSoR== 1 && guess_SSoR(1)>0 && guess_SSoR(2)>0  
    ReS=(guess_SSoR(1)*guess_SSoR(3)-guess_SSoR(2)*guess_SSoR(4))/(guess_SSoR(1)-guess_SSoR(2));
    lS=(A1*u1-guess_SSoR(1)*guess_SSoR(3))/(A1-guess_SSoR(1));
    c22=sqrt(k2*(m*guess_SSoR(5).^m+n*guess_SSoR(5).^(-n)));
    % distances 

    x1=lS*t;
    xReS = ReS*t;
    x22=(guess_SSoR(6)+c22)*t;
    x21 = x22;
    x2=(u2+c2)*t;
  
        % Plot A against x
        X=linspace(x1-1,x1,20); a=linspace(A1,A1,20); plot(X, a, 'color', "#f97306",LineWidth=1);
        hold on
        X=linspace(x1,x1,20); a=linspace(A1,guess_SSoR(1),20);plot(X, a, 'b',LineWidth=1);
         hold on
        X=linspace(x1,xReS,20); a=linspace(guess_SSoR(1),guess_SSoR(1),20);plot(X, a, 'g',LineWidth=1);
        hold on
        X=linspace(xReS,xReS,20); a=linspace(guess_SSoR(1), guess_SSoR(2),20); plot(X, a, 'color', "#7e1e9c",LineWidth=1);
        hold on
        X=linspace(xReS,x21,20); a=linspace(guess_SSoR(2), guess_SSoR(2),20); plot(X, a, 'c',LineWidth=1);   
        hold on
        X=linspace(x21,x22,20); a=linspace(guess_SSoR(2),A21,20);plot(X, a, 'm:',LineWidth=1);
        hold on
        X=linspace(x22,x2,20); a=linspace(A22, A2,20); plot(X, a, 'r',LineWidth=1);
        hold on
        X=linspace(x2,x2+1,20); a=linspace(A2, A2,20); plot(X, a, 'color', "#6e750e",LineWidth=1);
        legend({'$A_1$','1-$S$','$A_{11}$','$S_{\neq 0}$','$A_{12}$',...
                '3-$o$','2-$R$', '$A_2$'},...
                'interpreter','latex','FontSize',12,'Location','northeastoutside','Orientation','vertical');
    
    elseif exitflag==1  && x(3)>0 && x(1)>0 
          c1=sqrt(m*A1.^m+n*A1.^(-n)); 
          c11=sqrt(m*x(1).^m+n*x(1).^(-n)); 
          c22=sqrt(k2.*(m*(x(3)/A2bar).^m+n*(x(3)/A2bar).^(-n)));
          c2=sqrt(k2.*(m*(A2/A2bar).^m+n*(A2/A2bar).^(-n)));
          lS=(A1*u1-x(1)*x(2))/(A1-x(1));
          rS=(A2*u2-x(3)*x(4))/(A2-x(3));
          xS1=lS*t; 
          xS2=rS*t;
          x1=(u1-c1)*t;     
          x11=(x(2)-c11)*t;  
          x22=(x(4)+c22)*t;
          x2=(u2+c2)*t;

     if x(1) <= A1 && x(3) <= A2 %RoR

                X=linspace(x1-1,x1,20); a=linspace(A1,A1,20);plot(X, a, 'color', "#f97306",LineWidth=1);
                hold on
                X=linspace(x1,x11,20); a=linspace(A1,x(1),20);plot(X, a, 'b',LineWidth=1);
                hold on
                X=linspace(x11,0,20); a=linspace(x(1),x(1),20);plot(X, a, 'g',LineWidth=1);
                hold on
                X=linspace(0,0,20); a=linspace(x(1),x(3),20); plot(X,a,'m:',LineWidth=1);
                hold on
                X=linspace(0,x22,20); a=linspace(x(3),x(3),20);plot(X, a, 'c',LineWidth=1);
                hold on
                X=linspace(x22,x2,20); a=linspace(x(3),A2,20);plot(X, a, 'r',LineWidth=1);
                hold on
                X=linspace(x2,x2+1,20); a=linspace(A2,A2,20);  plot(X, a, 'color', "#6e750e",LineWidth=1);
                 %title('$RoR$','interpreter','latex')
                legend({'$A_1$','1-$R$','$A_{11}$', '3-$o$','$A_{22}$','2-$R$','$A_2$'},...
                    'interpreter','latex','FontSize',12,'Location','northeastoutside','Orientation','vertical');

     elseif x(1) <= A1 && x(3) > A2 %RoS

                X=linspace(x1-1,x1,20); a=linspace(A1,A1,20); plot(X, a, 'color', "#f97306",LineWidth=1);
                hold on
                X=linspace(x1,x11,20); a=linspace(A1,x(1),20);plot(X, a, 'b',LineWidth=1);
                hold on
                X=linspace(x11,0,20); a=linspace(x(1),x(1),20);plot(X, a, 'g',LineWidth=1);
                hold on
                X=linspace(0,0,20); a=linspace(x(1),x(3),20); plot(X, a,'m:',LineWidth=1);
                hold on
                X=linspace(0,xS2,20); a=linspace(x(3),x(3),20);plot(X, a, 'c',LineWidth=1);
                hold on
                X=linspace(xS2,xS2,20); a=linspace(x(3),A2,20);plot(X, a, 'r',LineWidth=1);
                hold on
                X=linspace(xS2,x2+1,20); a=linspace(A2,A2,20);plot(X, a,'color', "#6e750e",LineWidth=1);
                 %title('$RoS$','interpreter','latex')
                 legend({'$A_1$','1-$R$','$A_{11}$', '3-$o$','$A_{22}$','2-$S$','$A_2$'},...
                     'interpreter','latex','FontSize',12,'Location','northeastoutside','Orientation','vertical');

      elseif x(1) > A1 && x(3) <= A2  %SoR

                X=linspace(x1-1,xS1,20); a=linspace(A1,A1,20);plot(X, a, 'color', "#f97306",LineWidth=1);
                hold on
                X=linspace(xS1,xS1,20); a=linspace(A1,x(1),20);plot(X, a, 'b',LineWidth=1);
                hold on
                X=linspace(xS1,0,20); a=linspace(x(1),x(1),20);plot(X, a, 'g',LineWidth=1);
                hold on
                X=linspace(0,0,20); a=linspace(x(1),x(3),20);plot(X,a,'m:',LineWidth=1);
                hold on
                X=linspace(0,x22,20); a=linspace(x(3),x(3),20);plot(X, a, 'c',LineWidth=1);
                hold on
                X=linspace(x22,x2,20); a=linspace(x(3),A2,20); plot(X, a, 'r',LineWidth=1);
                hold on
                X=linspace(x2,x2+1,20); a=linspace(A2,A2,20);  plot(X, a, 'color', "#6e750e",LineWidth=1);
                %title('$SoR$','interpreter','latex')
                legend({'$A_1$','1-$S$','$A_{11}$', '3-$o$','$A_{22}$','2-$R$','$A_2$'},...
                    'interpreter','latex','FontSize',12,'Location','northeastoutside','Orientation','vertical');
                %xlim([-1.2,2.6])

       elseif x(1) > A1 && x(3) > A2  %SoS

                X=linspace(x1-1,xS1,20); a=linspace(A1,A1,20);plot(X, a, 'color', "#f97306",LineWidth=1);
                hold on
                X=linspace(xS1,xS1,20); a=linspace(A1,x(1),20);plot(X, a, 'b',LineWidth=1);
                hold on
                X=linspace(xS1,0,20); a=linspace(x(1),x(1),20);plot(X, a, 'g',LineWidth=1);
                hold on
                X=linspace(0,0,20); a=linspace(x(1),x(3),20);plot(X,a,'m:',LineWidth=1);
                hold on
                X=linspace(0,xS2,20); a=linspace(x(3),x(3),20);plot(X, a, 'c',LineWidth=1);
                hold on
                X=linspace(xS2,xS2,20); a=linspace(x(3),A2,20);plot(X, a, 'r',LineWidth=1);
                hold on
                X=linspace(xS2,x2+1,20); a=linspace(A2,A2,20);plot(X, a, 'color', "#6e750e",LineWidth=1);
                 %title('$SoS$','interpreter','latex')
                 legend({'$A_1$','1-$S$','$A_{11}$', '3-$o$','$A_{22}$','2-$S$','$A_2$'},...
                     'interpreter','latex','FontSize',12,'Location','northeastoutside','Orientation','vertical');
        
      end
   end
   xlabel('$x$','interpreter','latex','FontSize',12)
    ylabel('$A$','interpreter','latex','FontSize',12)
   ylim([0.2,1])
   %xlim([-1.2, 2.8])
   xlim([x1-0.5, x2+0.5])
   
 end
 