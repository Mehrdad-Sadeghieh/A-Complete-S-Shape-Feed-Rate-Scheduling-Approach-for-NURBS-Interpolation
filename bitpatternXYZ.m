function [BP_value ,BP,lengh1,lengh2,Time1] = bitpatternXYZ (X,Y,Ts)
%% 2 axis polar interpolation



BLU=91*pi/2000; %% mm

%% Drawing CNC Coordinate Systems

x0=1330;
y0=827;

lengh1=sqrt (X.^2 + (y0-Y).^2);
lengh2=sqrt ((X-x0).^2 + (y0-Y).^2);
% theta1=lengh1*360/2000/BLU;
% theta1=theta1-theta1(1);
% theta2=lengh2*360/2000/BLU;
% theta2=theta2-theta2(1);

%% Cartesian Coordinate System

% lengh1=X;
% lengh2=Y;
% theta1=lengh1*360/2000/BLU;
% theta1=theta1-theta1(1);
% theta2=lengh2*360/2000/BLU;
% theta2=theta2-theta2(1);


Time1=0:Ts:(size(X,1)-1)*Ts;

save('Time1','Time1');




%% creat intejer number

NN=size(X,1);

llengh1=(lengh1)./BLU;
llengh2=(lengh2)./BLU;

for i=1:NN-1
dlengh1(i)=(llengh1(i+1) - llengh1(i));
dlengh2(i)=(llengh2(i+1) - llengh2(i));
    
dlengh1_abs(i)=abs(llengh1(i+1) - llengh1(i));
dlengh2_abs(i)=abs(llengh2(i+1) - llengh2(i));
end

dLengh1=floor(dlengh1_abs);
dLengh2=floor(dlengh2_abs);

counter1=0;
counter2=0;

for i=1:NN-1
    counter1=counter1+dlengh1_abs(i)-floor(dlengh1_abs(i));
    counter2=counter2+dlengh2_abs(i)-floor(dlengh2_abs(i));

    if counter1>=1
        dLengh1(i)=dLengh1(i)+1;
        counter1=counter1-1;
    end
    if counter2>=1
        dLengh2(i)=dLengh2(i)+1;
        counter2=counter2-1;
    end

end

%% Creating BitPattern matrix

BP=zeros( NN-1,4);
for i=1:NN-1
if dlengh1(i)>0
    Bp(i,1)=dLengh1(i);
    Bp(i,2)=0;
else
    Bp(i,1)=0;
    Bp(i,2)=dLengh1(i);
end
if dlengh2(i)>0
    Bp(i,3)=dLengh2(i);
    Bp(i,4)=0;
else
    Bp(i,3)=0;
    Bp(i,4)=dLengh2(i);
end
end

BP=abs(Bp);

%% Comparing  

% LENGH1(1)=0;
% LENGH2(1)=0;
% 
%     for i=2:size(dLengh2,2)+1
%         LENGH1(i)=LENGH1(i-1)+BP(i-1,1)-BP(i-1,2);
%         LENGH2(i)=LENGH2(i-1)+BP(i-1,3)-BP(i-1,4);
%     end
%     
% LLENGH1(1)=0;
% LLENGH2(1)=0;
% 
%     for i=2:size(dlengh2,2)+1
%         LLENGH1(i)=LLENGH1(i-1)+dlengh1(i-1);
%         LLENGH2(i)=LLENGH2(i-1)+dlengh2(i-1);
%     end
    
% figure(8)
% plot(LENGH1,LENGH2)
% xlabel('LENGH1')
% ylabel('LENGH2')
% hold on
% plot(LLENGH1,LLENGH2,'r')



%% removing zero pulses when 4 axis are zero

% counter5=0;
% for i=1:NN-1
%     if BP(i,1)==0 && BP(i,2)==0 && BP(i,3)==0 && BP(i,4)==0 
%         counter5=counter5;
%     else
%         counter5=counter5+1;
%         BP_total(counter5,:)=BP(i,:);
%     end
% end
% BP_total=[BP_total;0 0 0 0];
% 
% BP=BP_total;

%%  binary definition for PCI-1240 input pulses
f=1;
h=1;
s=floor(size(BP,1)/16)+15;
BP_zp=zeros(1,s);
BP_zm=zeros(1,s);
BP_up=zeros(1,s);
BP_um=zeros(1,s);


for i=1:size(BP,1)
    
    BP_zp(f)=BP(i,1)*2^(h-1)+BP_zp(f);
    BP_zm(f)=BP(i,2)*2^(h-1)+BP_zm(f);
    BP_up(f)=BP(i,3)*2^(h-1)+BP_up(f);
    BP_um(f)=BP(i,4)*2^(h-1)+BP_um(f);
    h=h+1;
    
    if rem(i,16)==0
        f=f+1;
        h=1;
    end
end



BP_value=transpose([BP_zp ;BP_zm; BP_up; BP_um]);

%% save pulses for register programming

% save('BP_totaljerk' , 'BP_total')
% save('BP_valuejerk' , 'BP_value')

end
