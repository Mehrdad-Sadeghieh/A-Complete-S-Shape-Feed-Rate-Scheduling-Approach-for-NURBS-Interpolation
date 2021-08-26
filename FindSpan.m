function [spannumber]=FindSpan(n,p,u,U)

if u==U(n+1)
    spannumber=n;
else
    low=p;
    high=n+1;
    mid=floor((low+high)/2);
    while u<U(mid) || u>=U(mid+1)
        if u<U(mid)
            high=mid;
        else
            low=mid;
        end
        mid=floor((low+high)/2);
    end
    spannumber=mid;
end

    
