function [exitflag,guess] = resonance_RoSS_sing_gen(u1,u2,k2,m,n,A2bar,A1,A2,guess,A12_guess)

    exitflag = 0;
    
    % Find A12
    [A12,~,~]=fsolve(@find_A12,A12_guess);
    c12=sqrt((m*A12.^m+n*A12.^(-n)));
    u12=lRout(A12);
    
    % Define functions to find_A12 and LRout 
    function u=lRout(A)
       f1=@(a) (m*a.^(m-2)+n*a.^(-(n+2))).^(1/2);
        int1=quad(f1,A1,A);
    u=u1-int1;
    end

    function v=find_A12(A)
    c12=sqrt(m*(A)^m+n*(A)^(-n));
    u12=lRout(A);
    v=u12-c12;
    end

    % Initial celerities
    c1=sqrt(m*(A1/A2bar).^m+n*(A1/A2bar).^(-n));
    c2=sqrt(k2*(m*(A2/A2bar).^m+n*(A2/A2bar).^(-n)));
    
 % Define the cases of the system of nonlinear equations and apply the initial guess
    % Case 1
    [S,~,exit1]=fsolve(@RoSS1,guess);
    A21=S(1);
    A22=S(2);

    c12=sqrt(m*A12.^m+n*A12.^(-n));
    u21=A12*u12/A21;
    c21=sqrt(k2*(m*(A21/A2bar).^m+n*(A21/A2bar).^(-n)));
    u31=u21+sqrt(k2*(m*(A22.^(m+1)-A21.^(m+1))/((m+1)*A2bar.^m)+...
        n*(A21.^(-n+1)-A22.^(-n+1))/((n-1)*A2bar.^(-n)))*(A22-A21)/(A22*A21));

    u32=u2+sqrt(k2*(m*(A2.^(m+1)-A22.^(m+1))/((m+1)*A2bar.^m)+...
        n*(A22.^(-n+1)-A2.^(-n+1))/((n-1)*A2bar.^(-n)))*(A2-A22)/(A2*A22));
    u22=u32;

    sd1=(A22*u22-A21*u21)/(A22-A21);
    
    sd2=(A22*u22^2-A21*u21^2+k2*m/(m+1)/A2bar^m*(A22^(m+1)-A21^(m+1))+...
        k2*n/(n-1)/A2bar^n*(A21^(-n+1)-A22^(-n+1)))/(A22*u22-A21*u21);
       
    c22=sqrt(k2*(m*(A22/A2bar).^m+n*(A22/A2bar).^(-n))); 

    se1=(A2*u2-A22*u22)/(A2-A22);
    se2=(A2*u2^2-A22*u22^2+k2*m/(m+1)/A2bar^m*(A2^(m+1)-A22^(m+1))+...
        k2*n/(n-1)/A2bar^n*(A22^(-n+1)-A2^(-n+1)))/(A2*u2-A22*u22);
       
    if exit1>0&& (u12-c12)>(u1-c1) && sd1>0 && se1>0 && se1>sd1 && abs(sd1-sd2)<1e-4 && abs(se1-se2)<1e-4 && A1>A12 && (u2+c2)<se1 && se1<(u22+c22) &&...
            u21-c21>0 && isreal(A21) && isreal(A22) && (u22-c22)<sd1 && sd1<(u21-c21)
            disp(['exit1: ' num2str(exit1)]);
            exitflag = 1;
            guess = [A21,A22,u22,A12,u12,u21];
    end

    % Case 2
    [S,~,exit2]=fsolve(@RoSS2,guess);
    A21=S(1);
    A22=S(2);

    c12=sqrt(m*A12.^m+n*A12.^(-n));
    u21=A12*u12/A21;
    c21=sqrt(k2*(m*(A21/A2bar).^m+n*(A21/A2bar).^(-n)));
    u31=u21+sqrt(k2*(m*(A22.^(m+1)-A21.^(m+1))/((m+1)*A2bar.^m)+...
        n*(A21.^(-n+1)-A22.^(-n+1))/((n-1)*A2bar.^(-n)))*(A22-A21)/(A22*A21));

    u32=u2-sqrt(k2*(m*(A2.^(m+1)-A22.^(m+1))/((m+1)*A2bar.^m)+...
        n*(A22.^(-n+1)-A2.^(-n+1))/((n-1)*A2bar.^(-n)))*(A2-A22)/(A2*A22));
    u22=u32;

    sd1=(A22*u22-A21*u21)/(A22-A21);
    
    sd2=(A22*u22^2-A21*u21^2+k2*m/(m+1)/A2bar^m*(A22^(m+1)-A21^(m+1))+...
        k2*n/(n-1)/A2bar^n*(A21^(-n+1)-A22^(-n+1)))/(A22*u22-A21*u21);
       
    c22=sqrt(k2*(m*(A22/A2bar).^m+n*(A22/A2bar).^(-n))); 

    se1=(A2*u2-A22*u22)/(A2-A22);
    se2=(A2*u2^2-A22*u22^2+k2*m/(m+1)/A2bar^m*(A2^(m+1)-A22^(m+1))+...
        k2*n/(n-1)/A2bar^n*(A22^(-n+1)-A2^(-n+1)))/(A2*u2-A22*u22);
    
    if exit2>0 && (u12-c12)>(u1-c1) && sd1>0 && se1>0 && se1>sd1 && abs(sd1-sd2)<1e-4 && abs(se1-se2)<1e-4 && A1>A12 && (u2+c2)<se1 && se1<(u22+c22) &&...
            u21-c21>0 && isreal(A21) && isreal(A22) && (u22-c22)<sd1 && sd1<(u21-c21)
            disp(['exit2: ' num2str(exit2)]);
            exitflag = 1;
            guess = [A21,A22,u22,A12,u12,u21];
    end

    % Case 3
    [S,~,exit3]=fsolve(@RoSS3,guess);
    A21=S(1);
    A22=S(2);

    c12=sqrt(m*A12.^m+n*A12.^(-n));
    u21=A12*u12/A21;
    c21=sqrt(k2*(m*(A21/A2bar).^m+n*(A21/A2bar).^(-n)));
    u31=u21-sqrt(k2*(m*(A22.^(m+1)-A21.^(m+1))/((m+1)*A2bar.^m)+...
        n*(A21.^(-n+1)-A22.^(-n+1))/((n-1)*A2bar.^(-n)))*(A22-A21)/(A22*A21));

    u32=u2+sqrt(k2*(m*(A2.^(m+1)-A22.^(m+1))/((m+1)*A2bar.^m)+...
        n*(A22.^(-n+1)-A2.^(-n+1))/((n-1)*A2bar.^(-n)))*(A2-A22)/(A2*A22));
    u22=u32;

    sd1=(A22*u22-A21*u21)/(A22-A21);
    
    sd2=(A22*u22^2-A21*u21^2+k2*m/(m+1)/A2bar^m*(A22^(m+1)-A21^(m+1))+...
        k2*n/(n-1)/A2bar^n*(A21^(-n+1)-A22^(-n+1)))/(A22*u22-A21*u21);
       
    c22=sqrt(k2*(m*(A22/A2bar).^m+n*(A22/A2bar).^(-n))); 

    se1=(A2*u2-A22*u22)/(A2-A22);
    se2=(A2*u2^2-A22*u22^2+k2*m/(m+1)/A2bar^m*(A2^(m+1)-A22^(m+1))+...
        k2*n/(n-1)/A2bar^n*(A22^(-n+1)-A2^(-n+1)))/(A2*u2-A22*u22);
   
    if exit3>0 && (u12-c12)>(u1-c1) && sd1>0 && se1>0 && se1>sd1 && abs(sd1-sd2)<1e-4 && abs(se1-se2)<1e-4 && A1>A12 && (u2+c2)<se1 && se1<(u22+c22) &&...
            u21-c21>0 && isreal(A21) && isreal(A22) && (u22-c22)<sd1 && sd1<(u21-c21)
             disp(['exit3: ' num2str(exit3)]);
             exitflag = 1;
            guess = [A21,A22,u22,A12,u12,u21];
    end

    % Case 4
    [S,~,exit4]=fsolve(@RoSS4,guess);
    A21=S(1);
    A22=S(2);

    c12=sqrt(m*A12.^m+n*A12.^(-n));
    u21=A12*u12/A21;
    c21=sqrt(k2*(m*(A21/A2bar).^m+n*(A21/A2bar).^(-n)));
    u31=u21-sqrt(k2*(m*(A22.^(m+1)-A21.^(m+1))/((m+1)*A2bar.^m)+...
        n*(A21.^(-n+1)-A22.^(-n+1))/((n-1)*A2bar.^(-n)))*(A22-A21)/(A22*A21));

    u32=u2-sqrt(k2*(m*(A2.^(m+1)-A22.^(m+1))/((m+1)*A2bar.^m)+...
        n*(A22.^(-n+1)-A2.^(-n+1))/((n-1)*A2bar.^(-n)))*(A2-A22)/(A2*A22));
    u22=u32;

    sd1=(A22*u22-A21*u21)/(A22-A21);
    
    sd2=(A22*u22^2-A21*u21^2+k2*m/(m+1)/A2bar^m*(A22^(m+1)-A21^(m+1))+...
        k2*n/(n-1)/A2bar^n*(A21^(-n+1)-A22^(-n+1)))/(A22*u22-A21*u21);
       
    c22=sqrt(k2*(m*(A22/A2bar).^m+n*(A22/A2bar).^(-n))); 

    se1=(A2*u2-A22*u22)/(A2-A22);
    se2=(A2*u2^2-A22*u22^2+k2*m/(m+1)/A2bar^m*(A2^(m+1)-A22^(m+1))+...
        k2*n/(n-1)/A2bar^n*(A22^(-n+1)-A2^(-n+1)))/(A2*u2-A22*u22);
    
    if exit4>0 && (u12-c12)>(u1-c1) && sd1>0 && se1>0 && se1>sd1 && abs(sd1-sd2)<1e-4 && abs(se1-se2)<1e-4 && A1>A12 && (u2+c2)<se1 && se1<(u22+c22) &&...
            u21-c21>0 && isreal(A21) && isreal(A22) && (u22-c22)<sd1 && sd1<(u21-c21)
             disp(['exit4: ' num2str(exit4)]);
             exitflag = 1;
            guess = [A21,A22,u22,A12,u12,u21];
    end
    
    % Equations to solve
    function v=RoSS1(S)
        A21=S(1);
        A22=S(2);

        u21=A12*u12/A21;
        u31=u21+sqrt(k2*((m*(A22.^(m+1)-A21.^(m+1))/(m+1))+...
            n*(A21.^(-n+1)-A22.^(-n+1))/(n-1))*(A22-A21)/(A22*A21));
        u32=u2+sqrt(k2*(m*(A2.^(m+1)-A22.^(m+1))/(m+1)+...
             n*(A22.^(-n+1)-A2.^(-n+1))/(n-1))*(A2-A22)/(A2*A22));
        u22=u32;

        v(1)=u31-u32;
        v(2)=1/2*u12^2 + A12.^m-A12.^(-n)-1/2*u21^2-k2*((A21).^m-(A21).^(-n));
    end

    function v=RoSS2(S)
        A21=S(1);
        A22=S(2);

         u21=A12*u12/A21;
        u31=u21+sqrt(k2*((m*(A22.^(m+1)-A21.^(m+1))/(m+1))+...
            n*(A21.^(-n+1)-A22.^(-n+1))/(n-1))*(A22-A21)/(A22*A21));
        u32=u2-sqrt(k2*(m*(A2.^(m+1)-A22.^(m+1))/(m+1)+...
             n*(A22.^(-n+1)-A2.^(-n+1))/(n-1))*(A2-A22)/(A2*A22));
        u22=u32;

        v(1)=u31-u32;
        v(2)=1/2*u12^2 + A12.^m-A12.^(-n)-1/2*u21^2-k2*((A21).^m-(A21).^(-n));
    end

    function v=RoSS3(S)
        A21=S(1);
        A22=S(2);

         u21=A12*u12/A21;
        u31=u21-sqrt(k2*((m*(A22.^(m+1)-A21.^(m+1))/(m+1))+...
            n*(A21.^(-n+1)-A22.^(-n+1))/(n-1))*(A22-A21)/(A22*A21));
        u32=u2+sqrt(k2*(m*(A2.^(m+1)-A22.^(m+1))/(m+1)+...
             n*(A22.^(-n+1)-A2.^(-n+1))/(n-1))*(A2-A22)/(A2*A22));
        u22=u32;

        v(1)=u31-u32;
        v(2)=1/2*u12^2 + A12.^m-A12.^(-n)-1/2*u21^2-k2*((A21).^m-(A21).^(-n));
    end

    function v=RoSS4(S)
        A21=S(1);
        A22=S(2);

        u21=A12*u12/A21;
        u31=u21-sqrt(k2*((m*(A22.^(m+1)-A21.^(m+1))/(m+1))+...
            n*(A21.^(-n+1)-A22.^(-n+1))/(n-1))*(A22-A21)/(A22*A21));
        u32=u2-sqrt(k2*(m*(A2.^(m+1)-A22.^(m+1))/(m+1)+...
             n*(A22.^(-n+1)-A2.^(-n+1))/(n-1))*(A2-A22)/(A2*A22));
        u22=u32;

        v(1)=u31-u32;
        v(2)=1/2*u12^2 + A12.^m-A12.^(-n)-1/2*u21^2-k2*((A21).^m-(A21).^(-n));
    end

% Display the result
    disp('Resonant soln:');
    disp(guess);
    disp('RoSS');
end