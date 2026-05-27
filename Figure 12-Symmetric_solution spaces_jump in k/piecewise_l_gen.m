% Define the piece-wise function for the left state

function fl = piecewise_l_gen(x,A1,m,n)

    if x(1) <= A1       % left rarefaction
        fun = @(A) sqrt(m*A.^(m-2)+n*A.^(-n-2));
        fl = integral(fun,A1,x(1),'ArrayValued',true);

    else                % left shock
        g1=sqrt((x(1)-A1)/(A1*x(1)));
        h1=sqrt((m/(m+1)*(x(1).^(m+1)-A1.^(m+1)))- (n/(n-1)*(x(1).^(-n+1)-A1.^(-n+1))));
        fl = g1*h1;
    end
end
