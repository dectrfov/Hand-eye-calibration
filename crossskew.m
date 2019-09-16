function y=crossskew(a)
%y=[0 -a(3) a(2);a(3) 0 -a(1); -a(2) a(1) 0];
y=[a(3,2);a(1,3);a(2,1)];
end