% Define the piece-wise function for the right star region

function fr = piecewise_r_gen(x,A2,k2,m,n,A2bar)

    if x(3) <= A2       % right rarefaction
        fun = @(A) sqrt(k2*(m*(A/A2bar).^(m-2)+n*(A/A2bar).^(-n-2)));
        fr = integral(fun,A2,x(3),'ArrayValued',true);

    else                % right shock
        g=sqrt((x(3)-A2)/(x(3)*A2));
        h=sqrt(k2*((m/((m+1)*A2bar.^m)*(x(3).^(m+1)-A2.^(m+1)))- ...
            (n/((n-1)*A2bar.^(-n))*(x(3).^(-n+1)-A2.^(-n+1)))));
        fr = g*h;
    end
end