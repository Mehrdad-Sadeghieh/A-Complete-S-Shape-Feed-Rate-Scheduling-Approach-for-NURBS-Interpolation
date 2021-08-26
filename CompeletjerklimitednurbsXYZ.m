%%  A complete S-shape feed rate scheduling approach for NURBS interpolator  "Xu Du, Jie Huang, Li-Min Zhun"

clear
clc
close all
format long

%% Input variable
tic
load('try111111.txt');
Input=try111111;

An=10;          %% mm/s^2
Jn=30;         %% mm/s^3
At=An;
Jt=Jn;
Vmax=20;        %% mm/s
Vstart=0;
Vend=0;
Ts=0.001;       %% Second
chorderror=5e-6;
zita=1e-6;
TOL=1e-3;
initialvalue=[Jn An Vmax];

n=size(Input,1)-1;
p=3;
U=(0:1/(n-2):1);
U=[0 0 0  U(1:size(U,2)-1) 1 1 1 1];
% U=[0 0 0 1/4 1/4 2/4 2/4 3/4 3/4 1 1 1];  %%Circle
% U=[0 0 0 1 1 1];  %%Line
m=size(U,2)-1;
n=n+1;
m=m+1;


controlpoint_c=Input;
controlpoint(:,1)=controlpoint_c(:,1);
controlpoint(:,2)=controlpoint_c(:,2);
controlpoint(:,3)=controlpoint_c(:,3);

weight=zeros(1,n)+1;
% weight=[1,sqrt(2)/2,1,sqrt(2)/2,1,sqrt(2)/2,1,sqrt(2)/2,1];  %%Cricle

Ph=zeros(size(controlpoint,1),4);

for i=1:size(controlpoint,1)
    Ph(i,:)=[controlpoint(i,:) 1];
end
P=controlpoint;
Pw1=zeros(size(weight,2),4);
Pw2=zeros(size(weight,2),3);
for i=1:size(weight,2)
    Pw1(i,:)=Ph(i,:)*weight(i);
    Pw2(i,:)=P(i,:)*weight(i);
end
%%  calculation of critical points

d=3;
counter1=0;
counter2=0;
uu=0:0.0001:1;
maxcurvature_u=[];

Kcr1=(8*chorderror)/( (Vmax*Ts)^2+4*chorderror^2);
Kcr2=An/(Vmax^2);
Kcr3=sqrt(Jn/(Vmax^3));
kcr=[Kcr1 Kcr2 Kcr3];
Kcr=min(kcr);

C=zeros(10001,4);
curvature=zeros(1,10001);
curvature2=zeros(1,10001);


for u=0:0.0001:1
    counter1=counter1+1;
    CK= RatCurveDerivs(d,n,p,U,Pw2,u,weight);
    C(counter1,:)= CurvePoint(n,p,U,Pw1,u);
    curvature(counter1)=( norm( cross( CK(2,1:3),CK(3,1:3)))) / ( ( norm( CK(2,1:3)))^3);
    curvature2(counter1)=curvature(counter1);
    if counter1>=3
        if curvature(counter1-1)>=curvature(counter1) && curvature(counter1-1)>=curvature(counter1-2) && curvature(counter1-1)>=Kcr
            counter2=counter2+1;
            maxcurvature_u(counter2)=uu(counter1-1);
        end
    end
end

maxcurvature_u=[0 maxcurvature_u 1];
maxcurvature=zeros(1,size(maxcurvature_u,2));

for i=1:size(maxcurvature_u,2)
    CK= RatCurveDerivs(d,n,p,U,Pw2,maxcurvature_u(i),weight);
    maxcurvature(i)=( norm( cross( CK(2,1:3),CK(3,1:3)))) / ( ( norm( CK(2,1:3)))^3);
end

%% Arclength calculation

u=[maxcurvature_u;maxcurvature_u];
i=1;
counter3=1;
counter4=2;
counter5=0;
Arclength=zeros(1,counter2+1);
arclength=zeros(1,1);
uu=[u(2,1),u(2,2);u(2,1),u(2,2)];
uold=u(2,1);
unew=u(2,2);
while true
    uold=uu(2,i);
    if uold==maxcurvature_u(counter4) && counter5==counter3      
        Arclength(counter4-1)=arclength;
        counter4=counter4+1;
        
        if counter2+3==counter4
            break
        end
        arclength=0;
        uu=[u(2,counter4-1),u(2,counter4);u(2,counter4-1),u(2,counter4)];
        i=1;
    end
    uold=uu(2,i);
    unew=uu(2,i+1);
    
    umid1=(uu(2,i)+uu(2,i+1))/2;
    umid2=(uu(2,i)+3*uu(2,i+1))/4;
    umid3=(3*uu(2,i)+uu(2,i+1))/4;
    CK1= RatCurveDerivs(d,n,p,U,Pw2,uold,weight);
    CK2= RatCurveDerivs(d,n,p,U,Pw2,umid1,weight);
    CK3= RatCurveDerivs(d,n,p,U,Pw2,unew,weight);
    CK4= RatCurveDerivs(d,n,p,U,Pw2,umid2,weight);
    CK5= RatCurveDerivs(d,n,p,U,Pw2,umid3,weight);
    s1=( (unew-uold)/6)*( norm(CK1(2,1:3))+4*norm(CK2(2,1:3))+norm(CK3(2,1:3)));
    s21=( (unew-uold)/12)*( norm(CK1(2,1:3))+4*norm(CK5(2,1:3))+norm(CK2(2,1:3)));
    s22=( (unew-uold)/12)*( norm(CK2(2,1:3))+4*norm(CK4(2,1:3))+norm(CK3(2,1:3)));
    s2=s21+s22;
    
    if (1/10)*abs(s2-s1)<zita
        arclength=arclength+s2;
        counter5=counter3;
        i=i+1;
    else
        counter3=counter3+1;
        uu(2,1:i)=uu(1,1:i);
        uu(2,i+1)=(uu(1,i)+uu(1,i+1))/2;
        uu(2,(i+2):size(uu(1,:),2)+1)=uu(1,(i+1):size(uu(1,:),2));
        uu(1,:)=uu(2,:);
    end
    
end

% Arclength(counter4-1)=sqrt(400^2+200^2);
save('ArclengthXYZ','Arclength');
%%  Backward scanning

N=counter4-2;
Vback(N+1)=Vend;
for i=N:-1:1
    if i==1
        Vback(i)=Vstart;
    else
        x0=Vback(i+1)+TOL/5;
        F0=x0^3+Vback(i+1)*(x0)^2-((Vback(i+1))^2)*x0-(Vback(i+1))^3-Jt*Arclength(i)^2;
        dF0=3*x0^2+2*Vback(i+1)*x0-((Vback(i+1))^2);
        while true
            x=x0-(F0/dF0);
            if abs(x-x0)<TOL
                break
            else
                x0=x;
                F0=x0^3+Vback(i+1)*(x0)^2-((Vback(i+1))^2)*x0-(Vback(i+1))^3-Jt*Arclength(i)^2;
                dF0=3*x0^2+2*Vback(i+1)*x0-((Vback(i+1))^2);
            end
        end
        Vback(i)=x;
        
        
        if Vback(i)>((At^2/Jt)+Vback(i+1))
            y0=Vback(i)+TOL/5;
            g0=Jt*(y0^2)+(At^2)*y0-Jt*(Vback(i+1))^2+(At^2)*(Vback(i+1))-2*At*Jt*Arclength(i);
            dg0=2*Jt*y0+At^2;
            while true
                y=y0-(g0/dg0);
                if abs(y-y0)<TOL
                    break
                else
                    y0=y;
                    g0=Jt*(y0^2)+(At^2)*y0-Jt*(Vback(i+1))^2+(At^2)*(Vback(i+1))-2*At*Jt*Arclength(i);
                    dg0=2*Jt*y0+At^2;
                end
            end
            Vback(i)=y;
        end
        ucritical=maxcurvature_u(i);
        CK= RatCurveDerivs(d,n,p,U,Pw2,ucritical,weight);
        KCR=maxcurvature(i);%( norm( cross( CK(2,1:3),CK(3,1:3)))) / ( ( norm( CK(2,1:3)))^3);
        v1=Vback(i);
        v2=(2/Ts)*sqrt( (2*chorderror/KCR)-chorderror^2);
        v3=sqrt(An/KCR);
        v4=(Jn/(KCR^2))^(1/3);
        v5=Vmax;
        v=[v1 v2 v3 v4 v5];
        Vback(i)=min(v);
    end
end

%%  Forward scanning
Vcr=zeros(1,N+1);
Vcr(1)=Vstart;
Vforward=zeros(1,N+1);
Vforward(1)=Vstart;
for i=1:N
    if i==N
        Vcr(i+1)=Vend;
    else
        x1=Vcr(i)+TOL/5;
        F1=(x1^3+Vcr(i)*(x1)^2-((Vcr(i))^2)*x1-(Vcr(i))^3-Jt*Arclength(i)^2);
        dF1=(3*x1^2+2*Vcr(i)*x1-Vcr(i)^2);
        while true
            x=x1-(F1/dF1);
            if abs(x-x1)<TOL
                break
            else
                x1=x;
                F1=(x1^3+Vcr(i)*(x1)^2-((Vcr(i))^2)*x1-(Vcr(i))^3-Jt*Arclength(i)^2);
                dF1=(3*x1^2+2*Vcr(i)*x1-Vcr(i)^2);
            end
        end
        Vforward(i+1)=(x);
        
        if Vforward(i+1)>((At^2/Jt)+Vcr(i))%Vforward(i))
            y1=Vcr(i+1)+TOL/5;
            g1=(Jt*(y1^2)+(At^2)*y1-Jt*(Vcr(i))^2+(At^2)*(Vcr(i))-2*At*Jt*Arclength(i));
            dg1=(2*Jt*y1+At^2);
            while true
                y=y1-(g1/dg1);
                if abs(y-y1)<TOL
                    break
                else
                    y1=y;
                    g1=(Jt*(y1^2)+(At^2)*y1-Jt*(Vcr(i))^2+(At^2)*(Vcr(i))-2*At*Jt*Arclength(i));
                    dg1=(2*Jt*y1+At^2);
                end
            end
            Vforward(i+1)=(y);
        end
        
        Vcr(i+1)=min(Vforward(i+1),Vback(i+1));
    end
end

%%  S-shape feed rate profile generation

ssa=0;
ssd=0;
sla=0;
sld=0;
sscri=zeros(1,N);
slcri=zeros(1,N);
T1=zeros(1,N);
T2=zeros(1,N);
T3=zeros(1,N);
T4=zeros(1,N);
T5=zeros(1,N);
T6=zeros(1,N);
T7=zeros(1,N);
vmax=zeros(1,N);


for i=1:N
    sscri(i)=0;
    if Vcr(i)<Vcr(i+1)
        if (Vcr(i+1)-Vcr(i))<=((At^2)/Jt)
            ssa=(Vcr(i+1)+Vcr(i))*sqrt((Vcr(i+1)-Vcr(i))/Jt);
        end
        if (Vcr(i+1)-Vcr(i))>((At^2)/Jt)
            ssa=(0.5)*(Vcr(i+1)+Vcr(i))*((At/Jt)+(Vcr(i+1)-Vcr(i))/At);
        end
        sscri(i)=ssa;
    end
    
    
    if Vcr(i)>Vcr(i+1)
        if (Vcr(i)-Vcr(i+1))<=((At^2)/Jt)
            ssd=(Vcr(i+1)+Vcr(i))*sqrt((Vcr(i)-Vcr(i+1))/Jt);
        end
        if (Vcr(i)-Vcr(i+1))>((At^2)/Jt)
            ssd=(0.5)*(Vcr(i+1)+Vcr(i))*((At/Jt)+(Vcr(i)-Vcr(i+1))/At);
        end
        sscri(i)=ssd;
    end
    
    
    if (Vmax-Vcr(i))<=((At^2)/Jt)
        sla=(Vmax+Vcr(i))*sqrt((Vmax-Vcr(i))/Jt);
    end
    if (Vmax-Vcr(i))>((At^2)/Jt)
        sla=(0.5)*(Vmax+Vcr(i))*((At/Jt)+(Vmax-Vcr(i))/At);
    end
    if (Vmax-Vcr(i+1))<=((At^2)/Jt)
        sld=(Vcr(i+1)+Vmax)*sqrt((Vmax-Vcr(i+1))/Jt);
    end
    if (Vmax-Vcr(i+1))>((At^2)/Jt)
        sld=(0.5)*(Vcr(i+1)+Vmax)*((At/Jt)+(Vmax-Vcr(i+1))/At);
    end
    slcri(i)=sla+sld;
    
    T1(i)=0;
    T2(i)=0;
    T3(i)=0;
    T4(i)=0;
    T5(i)=0;
    T6(i)=0;
    T7(i)=0;
    
    
    if abs(Arclength(i)-sscri(i))<TOL      %% short block
        if Vcr(i)<Vcr(i+1) && (Vcr(i+1)-Vcr(i))<=((At^2)/Jt)
            T1(i)=sqrt((Vcr(i+1)-Vcr(i))/Jt); %  At/Jt;
            T3(i)=sqrt((Vcr(i+1)-Vcr(i))/Jt); %At/Jt;
            vmax(i)=Vcr(i+1);
        end
        if Vcr(i)<Vcr(i+1) && (Vcr(i+1)-Vcr(i))>((At^2)/Jt)
            vmax(i)=Vcr(i+1);
            T1(i)=At/Jt;
            T2(i)=((vmax(i)-Vcr(i))/At)-(At/Jt);
            T3(i)=At/Jt;
        end
        if Vcr(i)>Vcr(i+1) && (Vcr(i)-Vcr(i+1))<=((At^2)/Jt)
            vmax(i)=Vcr(i);
            T5(i)=sqrt((Vcr(i)-Vcr(i+1))/Jt);  %At/Jt;
            T7(i)=sqrt((Vcr(i)-Vcr(i+1))/Jt);   %At/Jt;
        end
        if Vcr(i)>Vcr(i+1) && (Vcr(i)-Vcr(i+1))>((At^2)/Jt)
            vmax(i)=Vcr(i);
            T5(i)=At/Jt;
            T6(i)=((vmax(i)-Vcr(i+1))/At)-(At/Jt);
            T7(i)=At/Jt;
        end
        
    elseif Arclength(i)>slcri(i)               %% long block
        
        vmax(i)=Vmax;
        
        if  Vcr(i)<vmax(i) && (vmax(i)-Vcr(i))<=((At^2)/Jt)
            T1(i)=sqrt((vmax(i)-Vcr(i))/Jt);  %At/Jt;
            T3(i)=sqrt((vmax(i)-Vcr(i))/Jt);  %At/Jt;
        end
        if Vcr(i)<vmax(i) && (vmax(i)-Vcr(i))>((At^2)/Jt)
            T1(i)=At/Jt;
            T2(i)=((vmax(i)-Vcr(i))/At)-(At/Jt);
            T3(i)=At/Jt;
        end
        if vmax(i)>Vcr(i+1) && (vmax(i)-Vcr(i+1))<=((At^2)/Jt)
            T5(i)=sqrt((vmax(i)-Vcr(i+1))/Jt);   %At/Jt;
            T7(i)=sqrt((vmax(i)-Vcr(i+1))/Jt);   %At/Jt;
        end
        if vmax(i)>Vcr(i+1) && (vmax(i)-Vcr(i+1))>((At^2)/Jt)
            T5(i)=At/Jt;
            T6(i)=((vmax(i)-Vcr(i+1))/At)-(At/Jt);
            T7(i)=At/Jt;
        end
        T4(i)=(Arclength(i)-sla-sld)/vmax(i);
        
    else                                                     %% medium block
        
        v0=max(Vcr(i+1),Vcr(i));
        v1=Vmax;
        v=(v0+v1)/2;
        sa=0;
        sd=0;
        while true
            if (v-Vcr(i))<=((At^2)/Jt)
                sa=(v+Vcr(i))*sqrt((v-Vcr(i))/Jt);
            end
            if (v-Vcr(i))>((At^2)/Jt)
                sa=(0.5)*(v+Vcr(i))*((At/Jt)+(v-Vcr(i))/At);
            end
            if (v-Vcr(i+1))<=((At^2)/Jt)
                sd=(Vcr(i+1)+v)*sqrt((v-Vcr(i+1))/Jt);
            end
            if (v-Vcr(i+1))>((At^2)/Jt)
                sd=(0.5)*(Vcr(i+1)+v)*((At/Jt)+(v-Vcr(i+1))/At);
            end
            if abs(sa+sd-Arclength(i))<=TOL
                vmax(i)=v;
                break
            end
            if (sa+sd-Arclength(i))>TOL
                v1=v;
                v=(v0+v)/2;
            end
            if (sa+sd-Arclength(i))<(-TOL)
                v0=v;
                v=(v+v1)/2;
            end
        end
        
        
        if Vcr(i)<vmax(i) && (vmax(i)-Vcr(i))<=((At^2)/Jt)
            T1(i)=sqrt((vmax(i)-Vcr(i))/Jt);   %At/Jt;
            T3(i)=sqrt((vmax(i)-Vcr(i))/Jt);   %At/Jt;
        end
        if Vcr(i)<vmax(i) && (vmax(i)-Vcr(i))>((At^2)/Jt)
            T1(i)=At/Jt;
            T2(i)=((vmax(i)-Vcr(i))/At)-(At/Jt);
            T3(i)=At/Jt;
        end
        if vmax(i)>Vcr(i+1) && (vmax(i)-Vcr(i+1))<=((At^2)/Jt)
            T5(i)=sqrt((vmax(i)-Vcr(i+1))/Jt);   %At/Jt;
            T7(i)=sqrt((vmax(i)-Vcr(i+1))/Jt);    %At/Jt;
        end
        if vmax(i)>Vcr(i+1) && (vmax(i)-Vcr(i+1))>((At^2)/Jt)
            T5(i)=At/Jt;
            T6(i)=((vmax(i)-Vcr(i+1))/At)-(At/Jt);
            T7(i)=At/Jt;
        end
    end
end

%%  Variable-jerk compensation

Na1=zeros(1,N);
Na2=zeros(1,N);
Na3=zeros(1,N);
Nc=zeros(1,N);
Nd1=zeros(1,N);
Nd2=zeros(1,N);
Nd3=zeros(1,N);
Jta=zeros(1,N);
Jtd=zeros(1,N);
Atnew=zeros(1,N);
Dtnew=zeros(1,N);

% % for i=1:N
% %     Tt(i)=T1(i)+T2(i)+T3(i)+T4(i)+T5(i)+T6(i)+T7(i);
% %     Nt(i)=ceil(Tt(i)/Ts);
% %     Ttptime(i)=Nt(i)*Ts;
% %     T1(i)= (Ttptime(i)/ Tt(i))*T1(i);
% %     T2(i)= (Ttptime(i)/ Tt(i))*T2(i);
% %     T3(i)= (Ttptime(i)/ Tt(i))*T3(i);
% %     T4(i)= (Ttptime(i)/ Tt(i))*T4(i);
% %     T5(i)= (Ttptime(i)/ Tt(i))*T5(i);
% %     T6(i)= (Ttptime(i)/ Tt(i))*T6(i);
% %     T7(i)= (Ttptime(i)/ Tt(i))*T7(i);
% % end


for i=1:N
    Na1(i)=ceil(T1(i)/Ts);
    Na2(i)=ceil(T2(i)/Ts);
    Na3(i)=ceil(T3(i)/Ts);
    Nc(i)=ceil(T4(i)/Ts);
    Nd1(i)=ceil(T5(i)/Ts);
    Nd2(i)=ceil(T6(i)/Ts);
    Nd3(i)=ceil(T7(i)/Ts);
end


for i=1:N
    if Na1(i)==0 && Na3(i)==0 && T1(i)~=0 && T3(i)~=0
        T1(i)=Ts;
        T3(i)=Ts;
        Na1(i)=1;
        Na3(i)=1;
    end
    
    if Nd1(i)==0 && Nd3(i)==0 && T7(i)~=0 && T5(i)~=0
        T5(i)=Ts;
        T7(i)=Ts;
        Nd1(i)=1;
        Nd3(i)=1;
    end
    
    if Na2(i)==0 && T2(i)~=0
        T2(i)=Ts;
        Na2(i)=1;
    end
    
    if Nd2(i)==0 && T6(i)~=0
        T6(i)=Ts;
        Nd2(i)=1;
    end
    
    if Nc(i)==0 && T4(i)~=0
        T4(i)=Ts;
        Nc(i)=1;
    end
end

for i=1:N
    T1(i)=Na1(i)*Ts;
    T2(i)=Na2(i)*Ts;
    T3(i)=Na3(i)*Ts;
    T4(i)=Nc(i)*Ts;
    T5(i)=Nd1(i)*Ts;
    T6(i)=Nd2(i)*Ts;
    T7(i)=Nd3(i)*Ts;
end

for i=1:N
    if T1(i)~=0 && T2(i)~=0 && T3(i)~=0 && T4(i)~=0 && T5(i)~=0 && T6(i)~=0 && T7(i)~=0 %%case1
        Jta(i)=(vmax(i)-Vcr(i))/(0.5*T1(i)^2+T1(i)*T2(i)+T1(i)*T3(i)-0.5*T3(i)^2);
        Jtd(i)=(vmax(i)-Vcr(i+1))/(-0.5*T5(i)^2+T5(i)*T6(i)+T5(i)*T7(i)+0.5*T7(i)^2);
        Atnew(i)=Jta(i)*T1(i);
        Dtnew(i)=Jtd(i)*T5(i);
    elseif T1(i)~=0 && T2(i)==0 && T3(i)~=0 && T4(i)~=0 && T5(i)~=0 && T6(i)~=0 && T7(i)~=0%%case2
        Jta(i)=(vmax(i)-Vcr(i))/(0.5*T1(i)^2+T1(i)*T3(i)-0.5*T3(i)^2);
        Jtd(i)=(vmax(i)-Vcr(i+1))/(-0.5*T5(i)^2+T5(i)*T6(i)+T5(i)*T7(i)+0.5*T7(i)^2);
        Atnew(i)=Jta(i)*T1(i);
        Dtnew(i)=Jtd(i)*T5(i);
    elseif T1(i)~=0 && T2(i)==0 && T3(i)~=0 && T4(i)~=0 && T5(i)~=0 && T6(i)==0 && T7(i)~=0%%case3
        Jta(i)=(vmax(i)-Vcr(i))/(0.5*T1(i)^2+T1(i)*T3(i)-0.5*T3(i)^2);
        Jtd(i)=(vmax(i)-Vcr(i+1))/(-0.5*T5(i)^2+T5(i)*T7(i)+0.5*T7(i)^2);
        Atnew(i)=Jta(i)*T1(i);
        Dtnew(i)=Jtd(i)*T5(i);
    elseif T1(i)~=0 && T2(i)~=0 && T3(i)~=0 && T4(i)~=0 && T5(i)~=0 && T6(i)==0 && T7(i)~=0%%case4
        Jta(i)=(vmax(i)-Vcr(i))/(0.5*T1(i)^2+T1(i)*T2(i)+T1(i)*T3(i)-0.5*T3(i)^2);
        Jtd(i)=(vmax(i)-Vcr(i+1))/(-0.5*T5(i)^2+T5(i)*T7(i)+0.5*T7(i)^2);
        Atnew(i)=Jta(i)*T1(i);
        Dtnew(i)=Jtd(i)*T5(i);
    elseif T1(i)==0 && T2(i)==0 && T3(i)==0 && T4(i)~=0 && T5(i)~=0 && T6(i)~=0 && T7(i)~=0%case5
        Jta(i)=0;
        Jtd(i)=(vmax(i)-Vcr(i+1))/(-0.5*T5(i)^2+T5(i)*T6(i)+T5(i)*T7(i)+0.5*T7(i)^2);
        Atnew(i)=Jta(i)*T1(i);
        Dtnew(i)=Jtd(i)*T5(i);
    elseif T1(i)~=0 && T2(i)~=0 && T3(i)~=0 && T4(i)~=0 && T5(i)==0 && T6(i)==0 && T7(i)==0%case6
        Jta(i)=(vmax(i)-Vcr(i))/(0.5*T1(i)^2+T1(i)*T2(i)+T1(i)*T3(i)-0.5*T3(i)^2);
        Jtd(i)=0;
        Atnew(i)=Jta(i)*T1(i);
        Dtnew(i)=Jtd(i)*T5(i);
    elseif T1(i)~=0 && T2(i)==0 && T3(i)~=0 && T4(i)~=0 && T5(i)==0 && T6(i)==0 && T7(i)==0%case7
        Jta(i)=(vmax(i)-Vcr(i))/(0.5*T1(i)^2+T1(i)*T3(i)-0.5*T3(i)^2);
        Jtd(i)=0;
        Atnew(i)=Jta(i)*T1(i);
        Dtnew(i)=Jtd(i)*T5(i);
    elseif T1(i)==0 && T2(i)==0 && T3(i)==0 && T4(i)~=0 && T5(i)~=0 && T6(i)==0 && T7(i)~=0%case8
        Jta(i)=0;
        Jtd(i)=(vmax(i)-Vcr(i+1))/(-0.5*T5(i)^2+T5(i)*T7(i)+0.5*T7(i)^2);
        Atnew(i)=Jta(i)*T1(i);
        Dtnew(i)=Jtd(i)*T5(i);
    elseif T1(i)==0 && T2(i)==0 && T3(i)==0 && T4(i)~=0 && T5(i)==0 && T6(i)==0 && T7(i)==0%case9
        Jta(i)=0;
        Jtd(i)=0;
        Atnew(i)=Jta(i)*T1(i);
        Dtnew(i)=Jtd(i)*T5(i);
    elseif T1(i)~=0 && T2(i)==0 && T3(i)~=0 && T4(i)==0 && T5(i)==0 && T6(i)==0 && T7(i)==0%case10
        Jta(i)=(Vcr(i+1)-Vcr(i))/(0.5*T1(i)^2+T1(i)*T3(i)-0.5*T3(i)^2);
        Jtd(i)=0;
        Atnew(i)=Jta(i)*T1(i);
        Dtnew(i)=Jtd(i)*T5(i);
    elseif T1(i)~=0 && T2(i)~=0 && T3(i)~=0 && T4(i)==0 && T5(i)==0 && T6(i)==0 && T7(i)==0%case11
        Jta(i)=(Vcr(i+1)-Vcr(i))/(0.5*T1(i)^2+T1(i)*T2(i)+T1(i)*T3(i)-0.5*T3(i)^2);
        Jtd(i)=0;
        Atnew(i)=Jta(i)*T1(i);
        Dtnew(i)=Jtd(i)*T5(i);
    elseif T1(i)==0 && T2(i)==0 && T3(i)==0 && T4(i)==0 && T5(i)~=0 && T6(i)==0 && T7(i)~=0%case12
        Jta(i)=0;
        Jtd(i)=(Vcr(i)-Vcr(i+1))/(-0.5*T5(i)^2+T5(i)*T7(i)+0.5*T7(i)^2);
        Atnew(i)=Jta(i)*T1(i);
        Dtnew(i)=Jtd(i)*T5(i);
    elseif T1(i)==0 && T2(i)==0 && T3(i)==0 && T4(i)==0 && T5(i)~=0 && T6(i)~=0 && T7(i)~=0%case13
        Jta(i)=0;
        Jtd(i)=(Vcr(i)-Vcr(i+1))/(-0.5*T5(i)^2+T5(i)*T6(i)+T5(i)*T7(i)+0.5*T7(i)^2);
        Atnew(i)=Jta(i)*T1(i);
        Dtnew(i)=Jtd(i)*T5(i);
    elseif T1(i)~=0 && T2(i)~=0 && T3(i)~=0 && T4(i)==0 && T5(i)~=0 && T6(i)~=0 && T7(i)~=0%case14
        Jta(i)=(vmax(i)-Vcr(i))/(0.5*T1(i)^2+T1(i)*T2(i)+T1(i)*T3(i)-0.5*T3(i)^2);
        Jtd(i)=(vmax(i)-Vcr(i+1))/(-0.5*T5(i)^2+T5(i)*T6(i)+T5(i)*T7(i)+0.5*T7(i)^2);
        Atnew(i)=Jta(i)*T1(i);
        Dtnew(i)=Jtd(i)*T5(i);
    elseif T1(i)~=0 && T2(i)~=0 && T3(i)~=0 && T4(i)==0 && T5(i)~=0 && T6(i)==0 && T7(i)~=0%case15
        Jta(i)=(vmax(i)-Vcr(i))/(0.5*T1(i)^2+T1(i)*T2(i)+T1(i)*T3(i)-0.5*T3(i)^2);
        Jtd(i)=(vmax(i)-Vcr(i+1))/(-0.5*T5(i)^2+T5(i)*T7(i)+0.5*T7(i)^2);
        Atnew(i)=Jta(i)*T1(i);
        Dtnew(i)=Jtd(i)*T5(i);
    elseif T1(i)~=0 && T2(i)==0 && T3(i)~=0 && T4(i)==0 && T5(i)~=0 && T6(i)~=0 && T7(i)~=0%case16
        Jta(i)=(vmax(i)-Vcr(i))/(0.5*T1(i)^2+T1(i)*T3(i)-0.5*T3(i)^2);
        Jtd(i)=(vmax(i)-Vcr(i+1))/(-0.5*T5(i)^2+T5(i)*T6(i)+T5(i)*T7(i)+0.5*T7(i)^2);
        Atnew(i)=Jta(i)*T1(i);
        Dtnew(i)=Jtd(i)*T5(i);
    elseif T1(i)~=0 && T2(i)==0 && T3(i)~=0 && T4(i)==0 && T5(i)~=0 && T6(i)==0 && T7(i)~=0%case17
        Jta(i)=(vmax(i)-Vcr(i))/(0.5*T1(i)^2+T1(i)*T3(i)-0.5*T3(i)^2);
        Jtd(i)=(vmax(i)-Vcr(i+1))/(-0.5*T5(i)^2+T5(i)*T7(i)+0.5*T7(i)^2);
        Atnew(i)=Jta(i)*T1(i);
        Dtnew(i)=Jtd(i)*T5(i);
    end
end
%%  Arc length error compensationNd1

e=zeros(1,N);
Ne1=zeros(1,N);
Ne2=zeros(1,N);
Ne3=zeros(1,N);
Ne5=zeros(1,N);
Ne6=zeros(1,N);
Ne7=zeros(1,N);
de_a=zeros(1,N);

for i=1:N
    vv1=Vcr(i)+(0.5)*Jta(i)*(Na1(i)*Ts)^2;
    vv2=vv1+Atnew(i)*(Na2(i)*Ts);
    %     vv3=vv2+Atnew(i)*(Na3(i)*Ts)-(0.5)*Jta(i)*(Na3(i)*Ts)^2;
    sa1=Vcr(i)*(Na1(i)*Ts)+(1/6)*Jta(i)*(Na1(i)*Ts)^3;
    sa2=sa1+vv1*(Na2(i)*Ts)+(0.5)*Atnew(i)*(Na2(i)*Ts)^2;
    saa=sa2+vv2*(Na3(i)*Ts)+(0.5)*Atnew(i)*(Na3(i)*Ts)^2-(1/6)*Jta(i)*(Na3(i)*Ts)^3;
    
    vv1=Vcr(i+1)+(0.5)*Jtd(i)*(Nd3(i)*Ts)^2;
    vv2=vv1+Dtnew(i)*(Nd2(i)*Ts);
    %     vv3=vv2+Dtnew(i)*(Nd1(i)*Ts)-(0.5)*Jtd(i)*(Nd1(i)*Ts)^2;
    sd1=Vcr(i+1)*(Nd3(i)*Ts)+(1/6)*Jtd(i)*(Nd3(i)*Ts)^3;
    sd2=sd1+vv1*(Nd2(i)*Ts)+(0.5)*Dtnew(i)*(Nd2(i)*Ts)^2;
    sdd=sd2+vv2*(Nd1(i)*Ts)+(0.5)*Dtnew(i)*(Nd1(i)*Ts)^2-(1/6)*Jtd(i)*(Nd1(i)*Ts)^3;
    
    sc=Nc(i)*Ts*vmax(i);
    
    e(i)=Arclength(i)-(saa+sdd+sc);
    
    Ne1(i)=(Na1(i)*(Na1(i)-1))/2;
    Ne2(i)=Na1(i)*Na2(i);
    Ne3(i)=(Na1(i)*(Na1(i)-1))/2;
    Ne5(i)=(Nd1(i)*(Nd1(i)-1))/2;
    Ne6(i)=Nd1(i)*Nd2(i);
    Ne7(i)=(Nd1(i)*(Nd1(i)-1))/2;
    
    if T1(i)~=0 && T2(i)~=0 && T3(i)~=0 && T4(i)~=0 && T5(i)~=0 && T6(i)~=0 && T7(i)~=0%case1
        de_a(i)=e(i)/(Ne1(i)+Ne2(i)+Ne3(i)+Ne5(i)+Ne6(i)+Ne7(i));
    elseif T1(i)~=0 && T2(i)==0 && T3(i)~=0 && T4(i)~=0 && T5(i)~=0 && T6(i)~=0 && T7(i)~=0%case2
        de_a(i)=e(i)/(Ne1(i)+Ne3(i)+Ne5(i)+Ne6(i)+Ne7(i));
    elseif T1(i)~=0 && T2(i)==0 && T3(i)~=0 && T4(i)~=0 && T5(i)~=0 && T6(i)==0 && T7(i)~=0%case3
        de_a(i)=e(i)/(Ne1(i)+Ne3(i)+Ne5(i)+Ne7(i));
    elseif T1(i)~=0 && T2(i)~=0 && T3(i)~=0 && T4(i)~=0 && T5(i)~=0 && T6(i)==0 && T7(i)~=0%case4
        de_a(i)=e(i)/(Ne1(i)+Ne2(i)+Ne3(i)+Ne5(i)+Ne7(i));
    elseif T1(i)==0 && T2(i)==0 && T3(i)==0 && T4(i)~=0 && T5(i)~=0 && T6(i)~=0 && T7(i)~=0%case5
        de_a(i)=e(i)/(Ne5(i)+Ne6(i)+Ne7(i));
    elseif T1(i)~=0 && T2(i)~=0 && T3(i)~=0 && T4(i)~=0 && T5(i)==0 && T6(i)==0 && T7(i)==0%case6
        de_a(i)=e(i)/(Ne1(i)+Ne2(i)+Ne3(i));
    elseif T1(i)~=0 && T2(i)==0 && T3(i)~=0 && T4(i)~=0 && T5(i)==0 && T6(i)==0 && T7(i)==0%case7
        de_a(i)=e(i)/(Ne1(i)+Ne3(i));
    elseif T1(i)==0 && T2(i)==0 && T3(i)==0 && T4(i)~=0 && T5(i)~=0 && T6(i)==0 && T7(i)~=0%case8
        de_a(i)=e(i)/(Ne5(i)+Ne7(i));
    elseif T1(i)==0 && T2(i)==0 && T3(i)==0 && T4(i)~=0 && T5(i)==0 && T6(i)==0 && T7(i)==0%case9
        de_a(i)=0;
    elseif T1(i)~=0 && T2(i)==0 && T3(i)~=0 && T4(i)==0 && T5(i)==0 && T6(i)==0 && T7(i)==0%case10
        de_a(i)=e(i)/(Ne1(i)+Ne3(i));
    elseif T1(i)~=0 && T2(i)~=0 && T3(i)~=0 && T4(i)==0 && T5(i)==0 && T6(i)==0 && T7(i)==0%case11
        de_a(i)=e(i)/(Ne1(i)+Ne2(i)+Ne3(i));
    elseif T1(i)==0 && T2(i)==0 && T3(i)==0 && T4(i)==0 && T5(i)~=0 && T6(i)==0 && T7(i)~=0%case12
        de_a(i)=e(i)/(Ne5(i)+Ne7(i));
    elseif T1(i)==0 && T2(i)==0 && T3(i)==0 && T4(i)==0 && T5(i)~=0 && T6(i)~=0 && T7(i)~=0%case13
        de_a(i)=e(i)/(Ne5(i)+Ne6(i)+Ne7(i));
    elseif T1(i)~=0 && T2(i)~=0 && T3(i)~=0 && T4(i)==0 && T5(i)~=0 && T6(i)~=0 && T7(i)~=0%case14
        de_a(i)=e(i)/(Ne1(i)+Ne2(i)+Ne3(i)+Ne5(i)+Ne6(i)+Ne7(i));
    elseif T1(i)~=0 && T2(i)~=0 && T3(i)~=0 && T4(i)==0 && T5(i)~=0 && T6(i)==0 && T7(i)~=0%case15
        de_a(i)=e(i)/(Ne1(i)+Ne2(i)+Ne3(i)+Ne5(i)+Ne7(i));
    elseif T1(i)~=0 && T2(i)==0 && T3(i)~=0 && T4(i)==0 && T5(i)~=0 && T6(i)~=0 && T7(i)~=0%case16
        de_a(i)=e(i)/(Ne1(i)+Ne3(i)+Ne5(i)+Ne6(i)+Ne7(i));
    elseif T1(i)~=0 && T2(i)==0 && T3(i)~=0 && T4(i)==0 && T5(i)~=0 && T6(i)==0 && T7(i)~=0%case17
        de_a(i)=e(i)/(Ne1(i)+Ne3(i)+Ne5(i)+Ne7(i));
    end
end

%%  Interpolation algorithm

counter5=2;
Ntt=zeros(1,N);

for i=1:N
    Ntt(i)=Na1(i)+Na2(i)+Na3(i)+Nc(i)+Nd1(i)+Nd2(i)+Nd3(i);
end

VVelocity=zeros(1,sum(Ntt)+1);
PPosition=zeros(1,sum(Ntt)+1);
VVelocity(1)=Vcr(1);
PPosition(1)=0;

for i=1:N
    constant=PPosition(counter5-1);
    for j=1:Na1(i)
        Jerk(counter5)=Jta(i);
        Acceleration(counter5)=Jta(i)*(j*Ts);
        VVelocity(counter5)=Vcr(i)+(0.5)*Jta(i)*(j*Ts)^2;
        PPosition(counter5)=constant+Vcr(i)*(j*Ts)+(1/6)*Jta(i)*(j*Ts)^3;
        counter5=counter5+1;
    end
    V1=VVelocity(counter5-1);
    Sa1=PPosition(counter5-1);
    
    for j=1:Na2(i)
        Jerk(counter5)=0;
        Acceleration(counter5)=Atnew(i);
        VVelocity(counter5)=V1+Atnew(i)*(j*Ts);
        PPosition(counter5)=Sa1+V1*(j*Ts)+(0.5)*Atnew(i)*(j*Ts)^2;
        counter5=counter5+1;
    end
    V2=VVelocity(counter5-1);
    Sa2=PPosition(counter5-1);
    
    for j=1:Na3(i)
        Jerk(counter5)=-Jta(i);
        Acceleration(counter5)=Atnew(i)-Jta(i)*(j*Ts);
        VVelocity(counter5)=V2+Atnew(i)*(j*Ts)-(0.5)*Jta(i)*(j*Ts)^2;
        PPosition(counter5)=Sa2+V2*(j*Ts)+(0.5)*Atnew(i)*(j*Ts)^2-(1/6)*Jta(i)*(j*Ts)^3;
        counter5=counter5+1;
    end
    V3=VVelocity(counter5-1);
    Sa=PPosition(counter5-1);
    
    for j=1:Nc(i)
        Jerk(counter5)=0;
        Acceleration(counter5)=0;
        VVelocity(counter5)=V3;
        PPosition(counter5)=Sa+V3*(j*Ts);
        counter5=counter5+1;
    end
    V4=VVelocity(counter5-1);
    Sac=PPosition(counter5-1);
    
    for j=1:Nd1(i)
        Jerk(counter5)=-Jtd(i);
        Acceleration(counter5)=-Jtd(i)*(j*Ts);
        VVelocity(counter5)=V4-(0.5)*Jtd(i)*(j*Ts)^2;
        PPosition(counter5)=Sac+V4*(j*Ts)-(1/6)*Jtd(i)*(j*Ts)^3;
        counter5=counter5+1;
    end
    V5=VVelocity(counter5-1);
    Sacd1=PPosition(counter5-1);
    
    for j=1:Nd2(i)
        Jerk(counter5)=0;
        Acceleration(counter5)=-Dtnew(i);
        VVelocity(counter5)=V5-Dtnew(i)*(j*Ts);
        PPosition(counter5)=Sacd1+V5*(j*Ts)-(0.5)*Dtnew(i)*(j*Ts)^2;
        counter5=counter5+1;
    end
    V6=VVelocity(counter5-1);
    Sacd2=PPosition(counter5-1);
    
    for j=1:Nd3(i)
        Jerk(counter5)=Jtd(i);
        Acceleration(counter5)=-Dtnew(i)+Jtd(i)*(j*Ts);
        VVelocity(counter5)=V6-Dtnew(i)*(j*Ts)+(0.5)*Jtd(i)*(j*Ts)^2;
        PPosition(counter5)=Sacd2+V6*(j*Ts)-(0.5)*Dtnew(i)*(j*Ts)^2+(1/6)*Jtd(i)*(j*Ts)^3;
        counter5=counter5+1;
    end
end
Jerk(end)=0;
dS=zeros(1,counter5-2);

for j=1:counter5-2
    dS(j)=PPosition(j+1)-PPosition(j);
    %     dV(j)=VVelocity(j+1)-VVelocity(j);
end

counter6=1;
dspa=zeros(1,sum(Ntt));
error=zeros(1,sum(Ntt));
hh=zeros(1,sum(Ntt));

for i=1:N
    for j=1:(Ntt(i))
        if j<=Na1(i)
            dspa(counter6)=dS(counter6)+(j-1)*de_a(i);
            error(counter6)=(j-1)*de_a(i);
            hh(counter6)=j-1;
            counter6=counter6+1;
        end
        if Na1(i)+1<=j && j<=Na1(i)+Na2(i)
            dspa(counter6)=dS(counter6)+Na1(i)*de_a(i);
            error(counter6)=Na1(i)*de_a(i);
            hh(counter6)=Na1(i);
            counter6=counter6+1;
        end
        if Na1(i)+Na2(i)+1<=j && j<=Na1(i)+Na2(i)+Na3(i)
            dspa(counter6)=dS(counter6)+(Na1(i)+Na2(i)+Na3(i)-j)*de_a(i);
            error(counter6)=(Na1(i)+Na2(i)+Na3(i)-j)*de_a(i);
            hh(counter6)=(Na1(i)+Na2(i)+Na3(i)-j);
            counter6=counter6+1;
        end
        if Na1(i)+Na2(i)+Na3(i)+1<=j && j<=Na1(i)+Na2(i)+Na3(i)+Nc(i)
            dspa(counter6)=dS(counter6);
            error(counter6)=0;
            hh(counter6)=0;
            counter6=counter6+1;
        end
        
        if Na1(i)+Na2(i)+Na3(i)+Nc(i)+1<=j && j<=Na1(i)+Na2(i)+Na3(i)+Nc(i)+Nd1(i)
            dspa(counter6)=dS(counter6)+(j-1-(Na1(i)+Na2(i)+Na3(i)+Nc(i)))*de_a(i);
            error(counter6)=(j-1-(Na1(i)+Na2(i)+Na3(i)+Nc(i)))*de_a(i);
            hh(counter6)=(j-1-(Na1(i)+Na2(i)+Na3(i)+Nc(i)));
            counter6=counter6+1;
        end
        
        if Na1(i)+Na2(i)+Na3(i)+Nc(i)+Nd1(i)+1<=j && j<=Na1(i)+Na2(i)+Na3(i)+Nc(i)+Nd1(i)+Nd2(i)
            dspa(counter6)=dS(counter6)+Nd1(i)*de_a(i);
            error(counter6)=Nd1(i)*de_a(i);
            hh(counter6)=Nd1(i);
            counter6=counter6+1;
        end
        
        if Na1(i)+Na2(i)+Na3(i)+Nc(i)+Nd1(i)+Nd2(i)+1<=j && j<=Na1(i)+Na2(i)+Na3(i)+Nc(i)+Nd1(i)+Nd2(i)+Nd3(i)
            dspa(counter6)=dS(counter6)+(Ntt(i)-j)*de_a(i);
            error(counter6)=(Ntt(i)-j)*de_a(i);
            hh(counter6)=(Ntt(i)-j);
            counter6=counter6+1;
        end
        
    end
end

Position=zeros(1,counter6);

for i=1:counter6-1
    Position(i+1)= Position(i)+dspa(i);
end
h=1;
h2=1;
sigmacount=zeros(1,counter6-1);
ua=zeros(1,counter6);
du=zeros(1,counter6-1);
ua1=0;
Ca=zeros(counter6,4);
Ca(1,:)= CurvePoint(n,p,U,Pw1,0);

d=3;
tic
for j=1:counter6-1
    
    CK2= RatCurveDerivs(d,n,p,U,Pw2,ua1,weight);
    
    sigma=norm(CK2(2,:));
    sigma1=dot(CK2(2,:),CK2(3,:))/sigma;

    ua2 = ua1+ (1/sigma) * dspa(j) - (0.5)*(sigma1/sigma^3)*dspa(j)^2;
    
    if ua2<0
        ua2=0;
        h2=h2+1;
    end
    if  ua2>1
        ua2=1;
        h=h+1;
        %         break
        
    end
    
    Ca(j+1,:)= CurvePoint(n,p,U,Pw1,ua2);
    ua1=ua2;
end
toc

TTP=0:Ts:(size(VVelocity,2)-1)*Ts;
TTE=0:Ts:(size(VVelocity,2)-2)*Ts;
TTEXYZ=TTE;
save('TTEXYZ','TTEXYZ');

[BP_value ,BP,length1,length2,Time1] = bitpatternXYZ (Ca(:,1),Ca(:,2),Ts);

%% Saving Outputs

save('Ferdowsi','length1')
save('Gear_length2','length2')
save('Gear_BP_value' , 'BP_value')
save('Gear_BP' , 'BP')
save('Gear_Ca' , 'Ca')

%% Figures

figure(1)
plot(C(:,1),C(:,2),'r','LineWidth',4)
hold on
plot(Ca(:,1),Ca(:,2),'k','LineWidth',4)
set(gca,'FontSize',30,'FontWeight','bold','FontName','Times New Roman')
xlabel('X','FontSize',30,'FontWeight','bold','FontName','Times New Roman')
ylabel('Y','FontSize',30,'FontWeight','bold','FontName','Times New Roman')


figure(2)
plot(TTP,VVelocity,'LineWidth',4)
set(gca,'FontSize',30,'FontWeight','bold','FontName','Times New Roman')
xlabel('Time','FontSize',30,'FontWeight','bold','FontName','Times New Roman')
ylabel('Velocity','FontSize',30,'FontWeight','bold','FontName','Times New Roman')

figure(3)
plot(TTE,error,'LineWidth',4)
xlabel('Time','FontSize',30,'FontWeight','bold','FontName','Times New Roman')
ylabel('Arclength Compensation Error','FontSize',30,'FontWeight','bold','FontName','Times New Roman')

figure(4)
plot(TTP,Position,'LineWidth',4)
xlabel('Time','FontSize',30,'FontWeight','bold','FontName','Times New Roman')
ylabel('Length','FontSize',30,'FontWeight','bold','FontName','Times New Roman')

figure(5)
plot(TTP,Acceleration,'LineWidth',4)
xlabel('Time','FontSize',30,'FontWeight','bold','FontName','Times New Roman')
ylabel('Acceleration','FontSize',30,'FontWeight','bold','FontName','Times New Roman')

figure(6)
plot(TTP,Jerk,'LineWidth',4)
xlabel('Time','FontSize',30,'FontWeight','bold','FontName','Times New Roman')
ylabel('Jerk','FontSize',30,'FontWeight','bold','FontName','Times New Roman')

figure(7)
plot(Ca(:,1),Ca(:,2),'k','LineWidth',4)
xlabel('X','FontSize',30,'FontWeight','bold','FontName','Times New Roman')
ylabel('Y','FontSize',30,'FontWeight','bold','FontName','Times New Roman')
hold on 
plot(Input(:,1),Input(:,2),'r--o','LineWidth',6)
legend('NURBS Curve','Control Points','FontSize',30,'FontWeight','bold','FontName','Times New Roman')

toc


