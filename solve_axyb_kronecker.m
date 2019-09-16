%A. Li, L. Wang, and D. Wu, "Simultaneous robot-world and hand-eye calibration using dual-quaternions and Kronecker product," Int. J. Phys. Sci., vol. 5, no. 10, pp. 1530¨C1536, 2010.
%input A (n*4)*4, B(n*4) *4
%output X 4*4 Y 4*4
function [X Y]=solve_axyb_kronecker(A,B)

AA=[];
YY=[];
for i=1:length(A)/4
RA=A(4*i-3:4*i,:);
RB=B(4*i-3:4*i,:);

AA=[AA;kron(RA(1:3,1:3),eye(3)) kron(-eye(3),RB(1:3,1:3)') zeros(9,6);zeros(3,9) kron(eye(3),RB(1:3,4)') -RA(1:3,1:3) eye(3)];
YY=[YY;zeros(9,1);RA(1:3,4)];
end

XX=inv(AA'*AA)*AA'*YY;
RX=[XX(1:3)';XX(4:6)';XX(7:9)'];
RY=[XX(10:12)';XX(13:15)';XX(16:18)'];
X=update_rot(RX);
X=[X XX(19:21);0 0 0 1];
Y=update_rot(RY);
Y=[Y XX(22:24);0 0 0 1];
end

function y=update_rot(x)
[aa bb cc]=svd(x);
if(det(aa*cc')<0)
    y=aa*[1 0 0;0 1 0; 0 0 -1]*cc';
else
    y=aa*cc';
end
end