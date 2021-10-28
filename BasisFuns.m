function [N] = BasisFuns(spannumber,u,p,U)



N(1)=1;
i=spannumber;
for j=1:p
    left(j) = u-U(i+1-j);
    right(j) = U(i+j)-u;
    saved = 0.0;
    r=1;
    while r<j+1
        temp = N(r)/(right(r)+left(j-r+1));
        N(r) = saved+right(r)*temp;
        saved = left(j-r+1)*temp;
        r=r+1;
    end
    N(j+1) = saved;
end
end