function [exitflag,guess] = resonance_SSoR_sing_gen(u1,u2,k2,m,n,A2bar,A1,A2,guess,A21_guess)
   exitflag = 0;
   options = optimoptions('fsolve', 'FunctionTolerance', 1e-8, 'OptimalityTolerance', 1e-8);

     % Find A21
    [A21,~,~]=fsolve(@find_A21,A21_guess,options);
    c21=sqrt(k2*(m*(A21/A2bar).^m+n*(A21/A2bar).^(-n)));
    u21=rRout(A21);

     % Define functions to find_A11 and rRout 
    function u=rRout(A)
       f1=@(a) (k2*(m*a.^(m-2)/A2bar.^m+n*a.^(-n-2)/A2bar.^(-n))).^(1/2);
      int1=integral(f1,A2,A);
    u=u2+int1;
    end

    function v=find_A21(A)
    c21=sqrt(k2*(m*(A/A2bar).^m+n*(A/A2bar).^(-n)));
    u21=rRout(A);
    v=u21+c21;
    end
     
   % Initial celerities
    c1=sqrt(m*A1.^m+n*A1.^(-n));
    c2=sqrt(k2.*(m*(A2/A2bar).^m+n*(A2/A2bar).^(-n)));

    % Case 1
    [S,~,exit1]=fsolve(@SSoR1,guess,options);
    A11=S(1);
    A12=S(2);

    u12=u21*A21/A12;
    c12=sqrt((m*A12.^m+n*A12.^(-n)));
    u31=u12+sqrt((m*(A11.^(m+1)-A12.^(m+1))/(m+1)+...
        n*(A12.^(-n+1)-A11.^(-n+1))/(n-1))*(A11-A12)/(A12*A11));

    u32=u1+sqrt((m*(A1.^(m+1)-A11.^(m+1))/(m+1)+...
        n*(A11.^(-n+1)-A1.^(-n+1))/(n-1))*(A1-A11)/(A1*A11));
    u11=u32;
    c11=sqrt((m*A11.^m+n*A11.^(-n)));
  
    sd1=(A11*u11-A12*u12)/(A11-A12);
    sd2=(A11*u11^2-A12*u12^2+ m/(m+1)*(A11^(m+1)-A12^(m+1))+...
        n/(n-1)*(A12^(-n+1)-A11^(-n+1)))/(A11*u11-A12*u12);

    se1=(A1*u1-A11*u11)/(A1-A11);
    se2=(A1*u1^2-A11*u11^2+m/(m+1)*(A1^(m+1)-A11^(m+1))+...
        n/(n-1)*(A11^(-n+1)-A1^(-n+1)))/(A1*u1-A11*u11);

    if exit1>0 && (u11-c11)<se1 && se1<(u1-c1)&& sd1<0 && se1<0 && se1<sd1 && abs(sd1-sd2)<1e-4 && abs(se1-se2)<1e-4 && abs(u31-u32)<1e-4 &&...
            (u12+c12)<0 && 0<(u2+c2) && isreal(A11) && isreal(A12)&& (u11+c11)>sd1 && sd1>(u12+c12)
        disp(['exit1: ' num2str(exit1)]); 
        exitflag = 1;
        guess = [A11,A12,u11,u12];
    end

    % Case 2
    [S,~,exit2]=fsolve(@SSoR2,guess,options);
    A11=S(1);
    A12=S(2);

    u12=u21*A21/A12;
    c12=sqrt((m*A12.^m+n*A12.^(-n)));

    u31=u12+sqrt((m*(A11.^(m+1)-A12.^(m+1))/(m+1)+...
        n*(A12.^(-n+1)-A11.^(-n+1))/(n-1))*(A11-A12)/(A12*A11));

    u32=u1-sqrt((m*(A1.^(m+1)-A11.^(m+1))/(m+1)+...
        n*(A11.^(-n+1)-A1.^(-n+1))/(n-1))*(A1-A11)/(A1*A11));
    u11=u32;
    c11=sqrt((m*A11.^m+n*A11.^(-n)));
  
    sd1=(A11*u11-A12*u12)/(A11-A12);
    sd2=(A11*u11^2-A12*u12^2+ m/(m+1)*(A11^(m+1)-A12^(m+1))+...
        n/(n-1)*(A12^(-n+1)-A11^(-n+1)))/(A11*u11-A12*u12);

    se1=(A1*u1-A11*u11)/(A1-A11);
    se2=(A1*u1^2-A11*u11^2+m/(m+1)*(A1^(m+1)-A11^(m+1))+...
        n/(n-1)*(A11^(-n+1)-A1^(-n+1)))/(A1*u1-A11*u11);
    
    if exit2>0 && (u11-c11)<se1 && se1<(u1-c1)&& sd1<0 && se1<0 && se1<sd1 && abs(sd1-sd2)<1e-4 && abs(se1-se2)<1e-4 && abs(u31-u32)<1e-4 &&...
            (u12+c12)<0 && 0<(u2+c2) && isreal(A11) && isreal(A12)&& (u11+c11)>sd1 && sd1>(u12+c12)
        disp(['exit2: ' num2str(exit2)]);
        exitflag = 1;
        guess = [A11,A12,u11,u12];
    end
    
    %case 3
    [S,~,exit3]=fsolve(@SSoR3,guess,options);
    A11=S(1);
    A12=S(2);

    u12=u21*A21/A12;
    c12=sqrt((m*A12.^m+n*A12.^(-n)));

    u31=u12-sqrt((m*(A11.^(m+1)-A12.^(m+1))/(m+1)+...
        n*(A12.^(-n+1)-A11.^(-n+1))/(n-1))*(A11-A12)/(A12*A11));

    u32=u1+sqrt((m*(A1.^(m+1)-A11.^(m+1))/(m+1)+...
        n*(A11.^(-n+1)-A1.^(-n+1))/(n-1))*(A1-A11)/(A1*A11));
    u11=u32;
    c11=sqrt((m*A11.^m+n*A11.^(-n)));
  
    sd1=(A11*u11-A12*u12)/(A11-A12);
    sd2=(A11*u11^2-A12*u12^2+ m/(m+1)*(A11^(m+1)-A12^(m+1))+...
        n/(n-1)*(A12^(-n+1)-A11^(-n+1)))/(A11*u11-A12*u12);

    se1=(A1*u1-A11*u11)/(A1-A11);
    se2=(A1*u1^2-A11*u11^2+m/(m+1)*(A1^(m+1)-A11^(m+1))+...
        n/(n-1)*(A11^(-n+1)-A1^(-n+1)))/(A1*u1-A11*u11);
   
    if exit3>0 && (u11-c11)<se1 && se1<(u1-c1)&& sd1<0 && se1<0 && se1<sd1 && abs(sd1-sd2)<1e-4 && abs(se1-se2)<1e-4 && abs(u31-u32)<1e-4 &&...
            (u12+c12)<0 && 0<(u2+c2) && isreal(A11) && isreal(A12)&& (u11+c11)>sd1 && sd1>(u12+c12)
        disp(['exit3: ' num2str(exit3)]);
        exitflag = 1;
        guess = [A11,A12,u11,u12];
    end

    % Case 4
    [S,~,exit4]=fsolve(@SSoR4,guess,options);
    A11=S(1);
    A12=S(2);

    u12=u21*A21/A12;
    c12=sqrt((m*A12.^m+n*A12.^(-n)));

    u31=u12-sqrt((m*(A11.^(m+1)-A12.^(m+1))/(m+1)+...
        n*(A12.^(-n+1)-A11.^(-n+1))/(n-1))*(A11-A12)/(A12*A11));

    u32=u1-sqrt((m*(A1.^(m+1)-A11.^(m+1))/(m+1)+...
        n*(A11.^(-n+1)-A1.^(-n+1))/(n-1))*(A1-A11)/(A1*A11));
    u11=u32;
    c11=sqrt((m*A11.^m+n*A11.^(-n)));
  
    sd1=(A11*u11-A12*u12)/(A11-A12);
    sd2=(A11*u11^2-A12*u12^2+ m/(m+1)*(A11^(m+1)-A12^(m+1))+...
        n/(n-1)*(A12^(-n+1)-A11^(-n+1)))/(A11*u11-A12*u12);

    se1=(A1*u1-A11*u11)/(A1-A11);
    se2=(A1*u1^2-A11*u11^2+m/(m+1)*(A1^(m+1)-A11^(m+1))+...
        n/(n-1)*(A11^(-n+1)-A1^(-n+1)))/(A1*u1-A11*u11);
    
    if exit4>0 && (u11-c11)<se1 && se1<(u1-c1)&& sd1<0 && se1<0 && se1<sd1 && abs(sd1-sd2)<1e-4 && abs(se1-se2)<1e-4 && abs(u31-u32)<1e-4 &&...
            (u12+c12)<0 && 0<(u2+c2) && isreal(A11) && isreal(A12)&& (u11+c11)>sd1 && sd1>(u12+c12)
        disp(['exit4: ' num2str(exit4)]);
        exitflag = 1;
        guess = [A11,A12,u11,u12];
    end
    
    % Equations to solve
    function v=SSoR1(S)
        A11=S(1);
        A12=S(2);

        u12=u21*A21/A12;
        u31=u12+sqrt((m*(A11.^(m+1)-A12.^(m+1))/(m+1)+...
            n*(A12.^(-n+1)-A11.^(-n+1))/(n-1))*(A11-A12)/(A12*A11));
        u32=u1+sqrt((m*(A1.^(m+1)-A11.^(m+1))/(m+1)+...
            n*(A11.^(-n+1)-A1.^(-n+1))/(n-1))*(A1-A11)/(A1*A11));

            v(1)=u31-u32;
            v(2)=1/2*u12^2 + A12.^m-A12.^(-n)-1/2*u21^2-k2*((A21/A2bar).^m-(A21/A2bar).^(-n));
    end

    function v=SSoR2(S)
        A11=S(1);
        A12=S(2);

        u12=u21*A21/A12;
        u31=u12+sqrt((m*(A11.^(m+1)-A12.^(m+1))/(m+1)+...
            n*(A12.^(-n+1)-A11.^(-n+1))/(n-1))*(A11-A12)/(A12*A11));
        u32=u1-sqrt((m*(A1.^(m+1)-A11.^(m+1))/(m+1)+...
            n*(A11.^(-n+1)-A1.^(-n+1))/(n-1))*(A1-A11)/(A1*A11));

            v(1)=u31-u32;
            v(2)=1/2*u12^2 + A12.^m-A12.^(-n)-1/2*u21^2-k2*((A21/A2bar).^m-(A21/A2bar).^(-n));
    end

    function v=SSoR3(S)
        A11=S(1);
        A12=S(2);

        u12=u21*A21/A12;
        u31=u12-sqrt((m*(A11.^(m+1)-A12.^(m+1))/(m+1)+...
            n*(A12.^(-n+1)-A11.^(-n+1))/(n-1))*(A11-A12)/(A12*A11));
        u32=u1+sqrt((m*(A1.^(m+1)-A11.^(m+1))/(m+1)+...
            n*(A11.^(-n+1)-A1.^(-n+1))/(n-1))*(A1-A11)/(A1*A11));

            v(1)=u31-u32;
            v(2)=1/2*u12^2 + A12.^m-A12.^(-n)-1/2*u21^2-k2*((A21/A2bar).^m-(A21/A2bar).^(-n));
    end

    function v=SSoR4(S)
        A11=S(1);
        A12=S(2);

        u12=u21*A21/A12;
        u31=u12-sqrt((m*(A11.^(m+1)-A12.^(m+1))/(m+1)+...
            n*(A12.^(-n+1)-A11.^(-n+1))/(n-1))*(A11-A12)/(A12*A11));
        u32=u1-sqrt((m*(A1.^(m+1)-A11.^(m+1))/(m+1)+...
            n*(A11.^(-n+1)-A1.^(-n+1))/(n-1))*(A1-A11)/(A1*A11));

            v(1)=u31-u32;
            v(2)=1/2*u12^2 + A12.^m-A12.^(-n)-1/2*u21^2-k2*((A21/A2bar).^m-(A21/A2bar).^(-n));
    end
    
% Display the result
    disp('Resonant soln:');
    disp(guess);
    disp('SSoR');
 end