function [exitflag,guess] = resonance_RoSR_sing_gen(u1,u2,k2,m,n,A2bar,A1,A2,guess, A12_guess)

    exitflag = 0;
 
      %  Find A12
    [A12,~,~]=fsolve(@find_A12,A12_guess);
    c12=sqrt((m*A12.^m+n*A12.^(-n)));
    u12=lRout(A12);

    % Define functions to find_A12 and lRout 
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

    % Case 1
    [S,~,exit1]=fsolve(@RoSR1,guess);
    A21=S(1);
    A22=S(2);

    c12=sqrt(m*A12.^m+n*A12.^(-n));
    u21=A12*u12/A21;
    c21=sqrt(k2*(m*(A21/A2bar).^m+n*(A21/A2bar).^(-n)));
    u31=u21+sqrt(k2*(m*(A22.^(m+1)-A21.^(m+1))/((m+1)*A2bar.^m)+...
        n*(A21.^(-n+1)-A22.^(-n+1))/((n-1)*A2bar.^(-n)))*(A22-A21)/(A22*A21));
    fun = @(A) sqrt(k2*(m*A.^(m-2)/A2bar.^(m)+n*A.^(-n-2)/A2bar.^(-n)));
       u32=u2+ quad(fun,A2,A22);
    u22=u32;
    c22=sqrt(k2*(m*(A22/A2bar).^m+n*(A22/A2bar).^(-n)));
    
    sd1=(A22*u22-A21*u21)/(A22-A21);
    sd2=(A22*u22^2-A21*u21^2+k2*m/(m+1)/A2bar^m*(A22^(m+1)-A21^(m+1))+...
        k2*n/(n-1)/A2bar^(-n)*(A21^(-n+1)-A22^(-n+1)))/(A22*u22-A21*u21);
    
    if exit1> 0 && (u12-c12)>(u1-c1) && sd1>0 && abs(sd1-sd2)<1e-4 && abs(u31-u32)<1e-4 && (u21-c21)>0 &&...
            (u22+c22)<(u2+c2) && isreal(A21) && isreal(A22)&& (u22-c22)<sd1 && sd1<(u21-c21)
        disp(['exit1: ' num2str(exit1)]);
        exitflag = 1;    
        guess = [A21,A22,A12];
    end

    % Case 2
    [S,~,exit2]=fsolve(@RoSR2,guess);
    A21=S(1);
    A22=S(2);
    
    c12=sqrt(m*A12.^m+n*A12.^(-n));
    u21=A12*u12/A21;
    c21=sqrt(k2*(m*(A21/A2bar).^m+n*(A21/A2bar).^(-n)));
    u31=u21-sqrt(k2*(m*(A22.^(m+1)-A21.^(m+1))/((m+1)*A2bar.^m)+...
        n*(A21.^(-n+1)-A22.^(-n+1))/((n-1)*A2bar.^(-n)))*(A22-A21)/(A22*A21));
    fun = @(A) sqrt(k2*(m*A.^(m-2)/A2bar.^(m)+n*A.^(-n-2)/A2bar.^(-n)));
       u32=u2+ quad(fun,A2,A22);
    u22=u32;
    c22=sqrt(k2*(m*(A22/A2bar).^m+n*(A22/A2bar).^(-n)));      
    
    sd1=(A22*u22-A21*u21)/(A22-A21);
    sd2=(A22*u22^2-A21*u21^2+k2*m/(m+1)/A2bar^m*(A22^(m+1)-A21^(m+1))+...
        k2*n/(n-1)/A2bar^(-n)*(A21^(-n+1)-A22^(-n+1)))/(A22*u22-A21*u21);

     if exit2>0 && (u12-c12)>(u1-c1) && sd1>0 && abs(sd1-sd2)<1e-4 && abs(u31-u32)<1e-4 && (u21-c21)>0 &&...
             (u22+c22)<(u2+c2) && isreal(A21) && isreal(A22)&& (u22-c22)<sd1 && sd1<(u21-c21)
         disp(['exit2: ' num2str(exit2)]);
         exitflag = 1;
         guess = [A21,A22,A12];   
     end
       

    % Equations to solve
    function v=RoSR1(S)
        A21=S(1);
        A22=S(2);
    
    u21=A12*u12/A21;
    u31=u21+sqrt(k2*(m*(A22.^(m+1)-A21.^(m+1))/((m+1)*A2bar.^m)+...
        n*(A21.^(-n+1)-A22.^(-n+1))/((n-1)*A2bar.^(-n)))*(A22-A21)/(A22*A21));
    fun = @(A) sqrt(k2*(m*A.^(m-2)/A2bar.^(m)+n*A.^(-n-2)/A2bar.^(-n)));
       u32=u2+ quad(fun,A2,A22);
    v(1)=u31-u32;        
    v(2)=1/2*u12^2 + A12.^m-A12.^(-n)-1/2*u21^2-k2*((A21/A2bar).^m-(A21/A2bar).^(-n));
    end
     
    function v=RoSR2(S)
        A21=S(1);
        A22=S(2);
   
    u21=A12*u12/A21;
    u31=u21-sqrt(k2*(m*(A22.^(m+1)-A21.^(m+1))/((m+1)*A2bar.^m)+...
        n*(A21.^(-n+1)-A22.^(-n+1))/((n-1)*A2bar.^(-n)))*(A22-A21)/(A22*A21));
    fun = @(A) sqrt(k2*(m*A.^(m-2)/A2bar.^(m)+n*A.^(-n-2)/A2bar.^(-n)));
       u32=u2+ quad(fun,A2,A22);
    v(1)=u31-u32;        
    v(2)=1/2*u12^2 + A12.^m-A12.^(-n)-1/2*u21^2-k2*((A21/A2bar).^m-(A21/A2bar).^(-n));
    end
    
 % Display the result
    disp('Resonant soln:');
    disp(guess);
    disp('RoSR');
end