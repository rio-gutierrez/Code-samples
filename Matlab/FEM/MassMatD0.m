function M = MassMatD0(x, p)
    %input mesh vector x and function p(x) to MassMatD0
    %output Mass Matrix M


    n = length(x)-1;         %number of subintervals
    M = zeros(n-1, n-1);     %allocate mass matrix
    %No need for half-hats due to vanishing BCs;
    %otherwise M would have dim (n+1)x(n+1)

    for i = 1:n-1
        h_minus = x(i+1) - x(i);                    %h_i       (index offset)
        xmid = (x(i+1) + x(i))/2;                   %m_i       (index offset)
        h_plus = x(i+2) - x(i+1);                   %h_{i+1}   (index offset)
        xmid_plus = (x(i+2) + x(i+1))/2;            %m_{i+1}   (index offset)

        M(i,i) = (h_minus/6) * ( p(xmid) + p(x(i+1)) )  + ...
                 (h_plus/6) * ( p(xmid_plus) + p(x(i+1)) ) ;

        if i ~= n-1
            M(i+1,i)  =  ( h_plus * p(xmid_plus) )/6;
            M(i,i+1) = M(i+1,i) ;
        end
    end

end
