%Polyfit function of ruler v. pixel data to second order polynomial
%excel file: Copy of Intial test 10_ruler variance - Repetative copies 3-26
%matlab file: polyfit.m - on desktop
%second degree polynomial, we get: y = ax^2 + bx + c
x = [30, 32.5, 35, 38, 47.5, 50, 57, 65]; %ruler position values
y = [75.025, 76.3, 78.2, 78.4, 79.925, 80.85, 79.533, 78.4]; % pixel values
p = polyfit(x,y,2);%2 means second order
f = polyval(p,x);
table = [x y f y-f];
plot(x,y,'o',x,f,'-')
axis([25 70 74 85])
x_new = (x - mean(x))/std(x); %This ensures that we are working with a small range of X_new values, [-1 1].
[p,S,mu] = polyfit(x,y,n);
A_new = p(1); B_new = p(2); C_new = p(3); mu1 = mu(1); mu2 = mu(2);
C=(A_new*(mu1)^2/(mu2)^2-(B_new*(mu1/mu2))+C_new);
B = (B_new*mu2 - 2*A_new*mu1)/(mu2)^2;
A = A_new/(mu2)^2;



