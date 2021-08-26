function [Nparametric] = BasisFunsparametric(spannumber,u,p,U)


i=spannumber;
for j=1:p
    left(j) = u-U(i+1-j);
    right(j) = U(i+j)-u;
    saved = 0.0;
    r=1;
    while r<j+1
        if j==1
            temp = 1/(right(r)+left(j-r+1));
        else
            temp = N(r)/(right(r)+left(j-r+1));
        end
        N(r)= saved+right(r)*temp;
        saved = left(j-r+1)*temp;
        r=r+1;
    end
    
    N(j+1)= saved;
end
Nparametric=sum(N);
end