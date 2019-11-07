%%
%Model-following with zero cancellation
clear ;
clc ;
format long;
m=100;%time(sec) ,sampling time=0.5(sec)
n=4;
theta(1,1)=0;%a1
theta(2,1)=0;%a2
theta(3,1)=0.01;%b0
theta(4,1)=0.2;%b1
phi(1,1)=0;
phi(1,2)=0;
phi(1,3)=0;
phi(1,4)=0;
w(1:50)=1;
w(51:100)=-1;
uc=[w w w w w];
uc=[uc uc];
p=[100 0 0 0 ;0 100 0 0; 0 0 1 0; 0 0 0 1];
lambda=0.8;
 %desired system:Am ,Bm
 am1=-1.3205;
 am2=0.4966;
 bm0=0.1761;
i=1;
 y(1)=0;
 u(1)=0;
for k=0:0.5:(m-0.5)     
    [p ,theta(:,i+1)]=rls_forgetting(p,theta(:,i),phi(i,:),y(i),lambda) ;
    %RLS:A,B
    a1=theta(1,i+1);
    a2=theta(2,i+1);
    b0=theta(3,i+1);
    b1=theta(4,i+1);
    %Controller:R,S,T
    r1=b1/b0;    
    s0=(am1-a1)/b0;
    s1=(am2-a2)/b0;
    t0=bm0/b0;    
    if i==1
        u(i+1)=-r1*u(i)+t0*uc(i)-s0*y(i)-s1*0;
        y(i+1)=1.6065*y(i)-0.6065*0+0.1065*u(i)+0.0902*0; %measured y
        phi(i+1,1)=-y(i);
        phi(i+1,2)=0;
        phi(i+1,3)=u(i);
        phi(i+1,4)=0;
    end
    if i==2
        u(i+1)=-r1*u(i)+t0*uc(i)-s0*y(i)-s1*y(i-1);
        y(i+1)=1.6065*y(i)-0.6065*y(i-1)+0.1065*u(i)+0.0902*u(i-1);      
        phi(i+1,1)=-y(i);
        phi(i+1,2)=-y(i-1);
        phi(i+1,3)=u(i);
        phi(i+1,4)=u(i-1);
    end    
    if i>=3               
        u(i+1)=-r1*u(i)+t0*uc(i)-s0*y(i)-s1*y(i-1);              
        y(i+1)=1.6065*y(i)-0.6065*y(i-1)+0.1065*u(i)+0.0902*u(i-1);     
        phi(i+1,1)=-y(i);
        phi(i+1,2)=-y(i-1);
        phi(i+1,3)=u(i);
        phi(i+1,4)=u(i-1);
    end           
    i=i+1;
end
figure(1)
subplot(211)
plot(0:0.5:m,y);
title('Model-following with zero cancellation')
xlabel('Time')
text(51, -0.08197,' y')
hold on
plot(0:0.5:m,uc(1:2*m+1));
text(21, -1,'uc')
axis([-inf, inf, -1.5, 1.5])
subplot(212)
plot(0:0.5:m,u);
xlabel('Time')
text(65, 1,' u')
axis([-inf, inf, -4.5, 4.5])

figure(2)
subplot(211)
plot(0:0.5:20,theta(1,1:41))
title('Model-following with zero cancellation')
xlabel('Time')
hold on
plot(0:0.5:20,theta(2,1:41))
hold on
plot([0,20],[-1.6065,-1.6065],'--')
hold on
plot([0,20],[0.6065,0.6065],'--')
text(3, -1.125,' a1')
text(3, -0.663,' a2')

subplot(212)
plot(0:0.5:20,theta(3,1:41))
xlabel('Time')
hold on
plot(0:0.5:20,theta(4,1:41))
hold on
plot([0,20],[0.1065,0.1065],'--')
hold on
plot([0,20],[0.0902,0.0902],'--')
text(2.5, 0.221,'b1')
text(1, 0.01,' b0')
%%
%Model-following without zero cancellation
clear ;
clc ;
format long;
m=100;%time(sec) ,sampling time=0.5(sec)
n=4;
theta(1,1)=0;%a1
theta(2,1)=0;%a2
theta(3,1)=0.01;%b0
theta(4,1)=0.2;%b1
phi(1,1)=0;
phi(1,2)=0;
phi(1,3)=0;
phi(1,4)=0;
w(1:50)=1;
w(51:100)=-1;
uc=[w w w w w];
uc=[uc uc];
p=[100 0 0 0 ;0 100 0 0; 0 0 1 0; 0 0 0 1];
lambda=0.8;
 %desired system:Am ,Bm
 am1=-1.3205;
 am2=0.4966;
 bm0=0.1761;   
i=1;
 y(1)=0;
 u(1)=1;
 a0=0;
for k=0:0.5:(m-0.5)     
    [p ,theta(:,i+1)]=rls_forgetting(p,theta(:,i),phi(i,:),y(i),lambda) ;

    %RLS: A,B
    a1=theta(1,i+1);
    a2=theta(2,i+1);
    b0=theta(3,i+1);
    b1=theta(4,i+1);
    %Controller:R,S,T
    
    r1=(a0*am2*b0*b0+(a2-am2-a0*am1)*b0*b1+(a0+am1-a1)*b1*b1)/(b1*b1-a1*b0*b1+a2*b0*b0);    
    s0=(b1*(a0*am1-a2-am1*a1+a1*a1+am2-a1*a0)+b0*(am1*a2-a1*a2-a0*am2+a0*a2))/(b1*b1-a1*b0*b1+a2*b0*b0);
    s1=(b1*(a1*a2-am1*a2+a0*am2-a0*a2)+b0*(a2*am2-a2*a2-a0*am2*a1+a0*a2*am1))/(b1*b1-a1*b0*b1+a2*b0*b0);
    t0=(1+am1+am2)/(b0+b1);        
    a0=(a2*r1+b1*s1)/am2;
    t1=t0*a0;
    if i==1
        u(i+1)=-r1*u(i)+t0*uc(i)+t1*0-s0*y(i)-s1*0;
        y(i+1)=1.6065*y(i)-0.6065*0+0.1065*u(i)+0.0902*0; %measured y
        phi(i+1,1)=-y(i);
        phi(i+1,2)=0;
        phi(i+1,3)=u(i);
        phi(i+1,4)=0;
    end
    if i==2
        u(i+1)=-r1*u(i)+t0*uc(i)+t1*uc(i-1)-s0*y(i)-s1*y(i-1);
        y(i+1)=1.6065*y(i)-0.6065*y(i-1)+0.1065*u(i)+0.0902*u(i-1);      
        phi(i+1,1)=-y(i);
        phi(i+1,2)=-y(i-1);
        phi(i+1,3)=u(i);
        phi(i+1,4)=u(i-1);
    end    
    if i>=3               
        u(i+1)=-r1*u(i)+t0*uc(i)+t1*uc(i-1)-s0*y(i)-s1*y(i-1);              
        y(i+1)=1.6065*y(i)-0.6065*y(i-1)+0.1065*u(i)+0.0902*u(i-1);     
        phi(i+1,1)=-y(i);
        phi(i+1,2)=-y(i-1);
        phi(i+1,3)=u(i);
        phi(i+1,4)=u(i-1);
    end           
    i=i+1;
end
figure(3)
subplot(211)
plot(0:0.5:m,y);
title('Model-following without zero cancellation')
xlabel('Time')
text(51, -0.08197,' y')
hold on
plot(0:0.5:m,uc(1:2*m+1));
text(22, -1,'uc')
axis([-inf, inf, -1.5, 1.5])
subplot(212)
plot(0:0.5:m,u);
xlabel('Time')
text(30, 1,' u')
axis([-inf, inf, -4.5, 4.5])

figure(4)
subplot(211)
plot(0:0.5:20,theta(1,1:41))
title('Model-following without zero cancellation')
xlabel('Time')
hold on
plot(0:0.5:20,theta(2,1:41))
hold on
plot([0,20],[-1.6065,-1.6065],'--')
hold on
plot([0,20],[0.6065,0.6065],'--')
text(3, -1.125,' a1')
text(3, 0,' a2')

subplot(212)
plot(0:0.5:20,theta(3,1:41))
xlabel('Time')
hold on
plot(0:0.5:20,theta(4,1:41))
hold on
plot([0,20],[0.1065,0.1065],'--')
hold on
plot([0,20],[0.0902,0.0902],'--')
text(2.5, 0.221,'b1')
text(1, 0.01,' b0')