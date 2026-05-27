% Define the piece-wise function for the right star region

function fr = piecewise_r_gen(x,A2,k2,m2,n2,A2bar)

    if x(3) <= A2       % right rarefaction
        fun = @(A) sqrt(k2*(m2*(A/A2bar).^(m2-2)+n2*(A/A2bar).^(-n2-2)));
        fr = integral(fun,A2,x(3),'ArrayValued',true);

    else                % right shock
        g=sqrt((x(3)-A2)/(x(3)*A2));
        h=sqrt(k2*((m2/((m2+1)*A2bar.^m2)*(x(3).^(m2+1)-A2.^(m2+1)))- ...
            (n2/((n2-1)*A2bar.^(-n2))*(x(3).^(-n2+1)-A2.^(-n2+1)))));
        fr = g*h;
    end
end