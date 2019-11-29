%Liao Wu and Hongliang Ren, Finding the Kinematic Base Frame of a Robot by Hand-Eye Calibration Using 3D Position Data, IEEE Transactions on Automation Science and Engineering, Vol. 14, No. 1, pp. 314-324, Jan. 2017.
%input A (n*4)*4, B(n*4) *4, X0,Y0 solved by closed form solution
%output X 4*4 Y 4*4
function [X,Y]=solve_axyb_guassnewton(X0,Y0,A,B)
tx=X0(1:3,4);
ry=Y0(1:3,1:3);
w0hat=logm(ry);
w0=crossskew(w0hat);
ty=Y0(1:3,4);

    for k=1:30
        Q=ry'*(eye(3)+((1-cos(norm(w0)))/norm(w0)^2)*w0hat+((norm(w0)-sin(norm(w0)))/norm(w0)^3)*w0hat^2); %eq.23
        J=[];
        f0=[];
        for i =1:length(A)/4
            RA=A(4*i-3:4*i,:);
            RB=B(4*i-3:4*i,:);
            J=[J;-ry*skewcross(RB(1:3,4))*Q ,eye(3) -RA(1:3,1:3)];
            f0=[f0;(ry*RB(1:3,4)+ty-RA(1:3,1:3)*tx-RA(1:3,4))]; %eq.23
        end
        y=-inv(J'*J)*J'*f0; %eq.27
        tx=tx+y(7:9);
        ty=ty+y(4:6);
        w0hat=w0hat+skewcross(y(1:3));
        w0=crossskew(w0hat);
        ry=expm(w0hat);
    end
Y(1:3,1:3)=ry;
X(1:3,1:3)=X0(1:3,1:3);

Y(1:3,4)=ty;
X(1:3,4)=tx;
X=[X;0 0 0 1];
Y=[Y;0 0 0 1];
end

function y=update_rot(x)
[aa bb cc]=svd(x);
if(det(aa*cc')<0)
    y=aa*[1 0 0;0 1 0; 0 0 -1]*cc';
else
    y=aa*cc';
end
end
