% Define the piece-wise function for the left state

function fl = piecewise_l_gen(x,A1,m1,n1)

    if x(1) <= A1       % left rarefaction
        fun = @(A) sqrt(m1*A.^(m1-2)+n1*A.^(-n1-2));
        fl = integral(fun,A1,x(1),'ArrayValued',true);

    else                % left shock
        g1=sqrt((x(1)-A1)/(A1*x(1)));
        h1=sqrt((m1/(m1+1)*(x(1).^(m1+1)-A1.^(m1+1)))- (n1/(n1-1)*(x(1).^(-n1+1)-A1.^(-n1+1))));
        fl = g1*h1;
    end
end
