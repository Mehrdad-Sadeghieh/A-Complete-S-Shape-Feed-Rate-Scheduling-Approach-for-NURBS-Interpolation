function nders = DersBasisFuns(spannumber,u,p,du,U)

ndu(1,1)=1;
i=spannumber;
for j=1:p
    left(j) = u-U(i+1-j);
    right(j) = U(i+j)-u;
    saved = 0;
    for r=1:j
        ndu(j+1,r)=right (r)+left (j-r+1) ;
        temp = ndu(r,j)/ndu(j+1,r);
        ndu(r,j+1) = saved+right(r)*temp;
        saved = left(j-r+1)*temp;
    end
    ndu(j+1,j+1) = saved;
end

for j=1:p+1
    nders(1,j)=ndu(j,p+1);
end




for r=1:p+1
    s1=1;
    s2=2;
    a(1,1)=1;
    for k=2:du+1
        d=0;
        rk=r-k;
        pk=p+1-k;
        if r>=k
            a(s2,1) = a(s1,1)/ndu(pk+1+1,rk+1);
            d = a(s2,1)*ndu(rk+1,pk+1);
        end
        if rk>=-1
            j1=1;
        else
            j1=-rk;
        end
        if (r-1-1)<=(pk)
            j2=k-1-1;
        else
            j2=p+1-r;
        end
        for j=j1:j2
            a(s2,j+1) = (a(s1,j+1)-a(s1,j))/ndu(pk+1+1,rk+j+1);
            d =d+ a(s2,j+1)*ndu(rk+j+1,pk+1);
        end
        if (r-1)<= (pk)
            a(s2,k) = -a(s1,k-1)/ndu(pk+1+1,r);
            d =d+ a(s2,k)*ndu(r,pk+1);
        end
        nders(k,r) =d;
        j=s1;
        s1=s2;
        s2=j;
    end
end
r=p;
for k=2:du+1
    for j=1:p+1
        nders(k,j) =nders(k,j)* r;
    end
    r =r* (p+1-k);
end

end

        
            
            
            
        
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
        
    






        