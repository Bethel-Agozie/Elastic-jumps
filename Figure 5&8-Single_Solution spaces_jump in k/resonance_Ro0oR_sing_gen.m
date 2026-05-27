function [exitflag,guess] = resonance_Ro0oR_sing_gen(u1,u2,k2,m,n,A2bar,A1,A2,guess,A11_guess)

    exitflag = 0;

    % Find A11
    [A11,~,~]=fsolve(@find_A11,A11_guess);
    c11=sqrt((m*A11.^m+n*A11.^(-n)));
    u11=lRout(A11);

     % Define functions to find_A11 and lRout 
    function u=lRout(A)
       f1=@(a) (m*a.^(m-2)+n*a.^(-(n+2))).^(1/2);
        int1=quad(f1,A1,A);
    u=u1-int1;
    end
  
    function v=find_A11(A)
    c11=sqrt(m*(A)^m+n*(A)^(-n));
    u11=lRout(A);
    v=u11-c11;
    end
    
    % Initial celerities
    c1=sqrt(m*A1.^m+n*A1.^(-n));
    c2=sqrt(k2.*(m*(A2/A2bar).^m+n*(A2/A2bar).^(-n)));

    % Define the cases of the system of nonlinear equations and apply the initial guess

    % Case 1
    [S,~,exit1]=fsolve(@Ro0oR1,guess);

    b=S(1);
    A12=S(2);
    A21=S(3);
    A22=S(4);
    
    kb=1+b*(k2-1); % kb=k1+b(k2-k1)

    u12=A11*u11/A12;
    c12=sqrt(kb*(m*A12.^m+n*A12.^(-n)));
    c11=sqrt(m*A11.^m+n*A11.^(-n));

    u21a=A12*u12/A21;
    u21b=u12+sqrt(kb*(m*(A21.^(m+1)-A12.^(m+1))/((m+1)*A2bar.^m)+...
        n*(A12.^(-n+1)-A21.^(-n+1))/((n-1)*A2bar.^(-n)))*(A21-A12)/(A21*A12));
    u21=u21a;
    c21=sqrt(kb*(m*(A21/A2bar).^m+n*(A21/A2bar).^(-n)));
    u22a=A21*u21/A22;
    fun = @(A) sqrt(k2*(m*A.^(m-2)/A2bar.^m+n*A.^(-n-2)/A2bar.^(-n)));
       u22b=u2+ quad(fun,A2,A22);
    u22=u22a;
    c22=sqrt(k2*(m*(A22/A2bar).^m+n*(A22/A2bar).^(-n)));

    if exit1 > 0 && (u11-c11)>(u1-c1) && u12-c12>0 && u21+c21>0 && A1>=A11 && A2>=A22 && (u22b+c22) < (u2+c2) &&...
            b<=1 && b>=0 && A11<=A22 && isreal(A12) && isreal(A21) && isreal(A22)&& isreal(b)   
        disp(['exit1: ' num2str(exit1)]);
        exitflag = 1;
        guess = [b,A12,A21,A22];
    end

    % Case2
    [S,~,exit2]=fsolve(@Ro0oR2,guess);
    b=S(1);
    A12=S(2);
    A21=S(3);
    A22=S(4);
    
    kb=1+b*(k2-1); % kb=k1+b(k2-k1)

    u12=A11*u11/A12;
    c12=sqrt(kb*(m*A12.^m+n*A12.^(-n)));
    c11=sqrt(m*A11.^m+n*A11.^(-n));

    u21a=A12*u12/A21;
    u21b=u12-sqrt(kb*(m*(A21.^(m+1)-A12.^(m+1))/((m+1)*A2bar.^m)+...
        n*(A12.^(-n+1)-A21.^(-n+1))/((n-1)*A2bar.^(-n)))*(A21-A12)/(A21*A12));
    u21=u21a;
    c21=sqrt(kb*(m*(A21/A2bar).^m+n*(A21/A2bar).^(-n)));
    u22a=A21*u21/A22;
    fun = @(A) sqrt(k2*(m*A.^(m-2)/A2bar.^m+n*A.^(-n-2)/A2bar.^(-n)));
       u22b=u2+ quad(fun,A2,A22);
    u22=u22a;
    c22=sqrt(k2*(m*(A22/A2bar).^m+n*(A22/A2bar).^(-n)));

    if exit2 > 0 && (u11-c11)>(u1-c1) && u12-c12>0 && u21+c21>0 && A1>=A11 && A2>=A22 && (u22b+c22) < (u2+c2) &&...
            b<=1 && b>=0 && A11<=A22 && isreal(A12) && isreal(A21) && isreal(A22)&& isreal(b)
       disp(['exit2: ' num2str(exit2)]);
       exitflag = 1;
       guess = [b,A12,A21,A22];
    end
   

   % Equations to solve
   function v=Ro0oR1(S)
    b=S(1);
    A12=S(2);
    A21=S(3);
    A22=S(4);

    kb=1+b*(k2-1); % kb=k1+b(k2-k1)

    u12=A11*u11/A12;
    u21a=A12*u12/A21;
    u21b=u12+sqrt(kb*(m*(A21.^(m+1)-A12.^(m+1))/((m+1)*A2bar.^m)+...
        n*(A12.^(-n+1)-A21.^(-n+1))/((n-1)*A2bar.^(-n)))*(A21-A12)/(A21*A12));
    u21=u21a;
    u22a=A21*u21/A22;
    fun = @(A) sqrt(k2*(m*A.^(m-2)/A2bar.^m+n*A.^(-n-2)/A2bar.^(-n)));
       u22b=u2+ quad(fun,A2,A22);
    u22=u22a;

        v(1)=1/2*u11^2+A11.^m-A11.^(-n)-1/2*u12^2-kb*(A12.^m-A12.^(-n));
        v(2)=u21a-u21b;
        v(3)=u22a-u22b;
        v(4)=1/2*u21^2+kb*((A21/A2bar).^m-(A21/A2bar).^(-n))-1/2*u22^2-k2*((A22/A2bar).^m-(A22/A2bar).^(-n));
    end

    function v=Ro0oR2(S)
    b=S(1);
    A12=S(2);
    A21=S(3);
    A22=S(4);
    
    kb=1+b*(k2-1); % kb=k1+b(k2-k1)

    u12=A11*u11/A12;
    u21a=A12*u12/A21;
    u21b=u12-sqrt(kb*(m*(A21.^(m+1)-A12.^(m+1))/((m+1)*A2bar.^m)+...
        n*(A12.^(-n+1)-A21.^(-n+1))/((n-1)*A2bar.^(-n)))*(A21-A12)/(A21*A12));
    u21=u21a;
    u22a=A21*u21/A22;
    fun = @(A) sqrt(k2*(m*A.^(m-2)/A2bar.^m+n*A.^(-n-2)/A2bar.^(-n)));
       u22b=u2+ quad(fun,A2,A22);
    u22=u22a;

        v(1)=1/2*u11^2+A11.^m-A11.^(-n)-1/2*u12^2-kb*(A12.^m-A12.^(-n));
        v(2)=u21a-u21b;
        v(3)=u22a-u22b;
        v(4)=1/2*u21^2+kb*((A21/A2bar).^m-(A21/A2bar).^(-n))-1/2*u22^2-k2*((A22/A2bar).^m-(A22/A2bar).^(-n));
    end

 % Display the result
    disp('Resonant soln:');
    disp(guess);
    disp('Ro0oR');
end
