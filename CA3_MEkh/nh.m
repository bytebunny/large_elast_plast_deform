%**************************************************************
%* AceGen    6.503 Windows (12 Sep 15)                        *
%*           Co. J. Korelc  2013           1 Mar 16 11:28:09  *
%**************************************************************
% User     : Full professional version
% Notebook : nh
% Evaluation time                 : 3 s     Mode  : Optimal
% Number of formulae              : 113     Method: Automatic
% Subroutine                      : nh size: 1735
% Total size of Mathematica  code : 1735 subexpressions
% Total size of Matlab code      : 5005 bytes 

%*********************** F U N C T I O N **************************
function[S,dSdE]=nh(CC,mu,lambda);
persistent v;
if size(v)<212
  v=zeros(212,'double');
end;
v(181)=CC(2)*CC(3)-CC(5)*CC(9);
v(180)=CC(4)*CC(6)-CC(1)*CC(9);
v(179)=-(CC(2)*CC(6))+CC(8)*CC(9);
v(201)=2e0*v(179);
v(178)=-(CC(3)*CC(4))+CC(7)*CC(9);
v(194)=mu*v(178);
v(191)=2e0*v(178);
v(177)=CC(1)*CC(2)-CC(4)*CC(8);
v(176)=CC(5)*CC(6)-CC(3)*CC(8);
v(175)=-(CC(1)*CC(5))+CC(7)*CC(8);
v(199)=mu*v(175);
v(196)=2e0*v(175);
v(174)=CC(1)*CC(3)-CC(6)*CC(7);
v(173)=CC(4)*CC(5)-CC(2)*CC(7);
v(14)=CC(6)*v(173)+CC(9)*v(175)+CC(3)*v(177);
v(182)=mu/v(14);
v(52)=log(sqrt(v(14)));
v(183)=lambda*v(52);
v(43)=1/Power(v(14),2);
v(51)=-(v(175)*v(43));
v(50)=-(v(178)*v(43));
v(49)=-(v(179)*v(43));
v(48)=-(v(173)*v(43));
v(192)=-(mu*v(48));
v(47)=-(v(180)*v(43));
v(190)=-(mu*v(47));
v(46)=-(v(176)*v(43));
v(204)=mu*v(46);
v(45)=-(v(177)*v(43));
v(44)=-(v(174)*v(43));
v(42)=-(v(181)*v(43));
v(122)=CC(7)*v(182);
v(117)=CC(8)*v(182);
v(114)=CC(9)*v(182);
v(104)=CC(4)*v(182);
v(96)=CC(6)*v(182);
v(92)=CC(5)*v(182);
v(54)=(lambda*v(43))/2e0;
v(61)=v(183)*v(50)+v(178)*v(54);
v(60)=v(183)*v(49)+v(179)*v(54);
v(59)=v(183)*v(48)+v(173)*v(54);
v(197)=v(192)+v(59);
v(58)=v(183)*v(47)+v(180)*v(54);
v(202)=v(190)+v(58);
v(57)=v(183)*v(46)+v(176)*v(54);
v(200)=-v(204)+v(57);
v(56)=v(183)*v(45)+v(177)*v(54);
v(188)=-(mu*v(45))+v(56);
v(189)=2e0*v(188);
v(55)=v(183)*v(44)+v(174)*v(54);
v(187)=-(mu*v(44))+v(55);
v(186)=2e0*v(187);
v(53)=v(183)*v(42)+v(181)*v(54);
v(184)=-(mu*v(42))+v(53);
v(185)=2e0*v(184);
v(17)=v(183)/v(14);
v(207)=v(17)-v(182);
v(123)=-(CC(7)*v(17));
v(195)=v(122)+v(123);
v(118)=-(CC(8)*v(17));
v(198)=v(117)+v(118);
v(115)=-(CC(9)*v(17));
v(193)=v(114)+v(115);
v(110)=CC(2)*v(17);
v(105)=-(CC(4)*v(17));
v(205)=v(104)+v(105);
v(100)=CC(1)*v(17);
v(97)=-(CC(6)*v(17));
v(206)=v(96)+v(97);
v(93)=-(CC(5)*v(17));
v(203)=v(92)+v(93);
v(87)=CC(3)*v(17);
v(88)=-(v(194)*v(46));
v(101)=v(175)*v(190);
v(111)=v(179)*v(192);
v(127)=2e0*(v(174)*v(53)+v(87)+v(88));
v(128)=2e0*(v(110)+v(111)+v(177)*v(53));
v(129)=v(176)*v(185);
v(130)=2e0*(v(180)*v(184)+v(193));
v(131)=v(173)*v(185);
v(132)=v(179)*v(185);
v(133)=v(178)*v(185);
v(134)=2e0*(v(175)*v(184)+v(203));
v(136)=2e0*(v(100)+v(101)+v(177)*v(55));
v(137)=v(176)*v(186);
v(138)=v(180)*v(186);
v(139)=2e0*(v(173)*v(187)+v(195));
v(140)=2e0*(v(179)*v(187)+v(206));
v(141)=v(178)*v(186);
v(142)=v(175)*v(186);
v(144)=2e0*(v(176)*v(188)+v(198));
v(145)=v(180)*v(189);
v(146)=v(173)*v(189);
v(147)=v(179)*v(189);
v(148)=2e0*(v(178)*v(188)+v(205));
v(149)=v(175)*v(189);
v(150)=2e0*(CC(3)*v(182)+v(178)*v(57)-v(87)+v(88));
v(151)=v(191)*v(202);
v(152)=v(191)*v(197);
v(153)=-2e0*(v(193)+v(194)*v(49)-v(178)*v(60));
v(155)=-2e0*(v(195)+v(199)*v(50)-v(175)*v(61));
v(156)=v(196)*v(200);
v(157)=2e0*(-v(100)+v(101)+CC(1)*v(182)+v(175)*v(58));
v(158)=v(196)*v(197);
v(159)=-2e0*(v(198)+v(199)*v(49)-v(175)*v(60));
v(161)=v(200)*v(201);
v(162)=v(201)*v(202);
v(163)=2e0*(-v(110)+v(111)+CC(2)*v(182)+v(179)*v(59));
v(165)=-2e0*(v(203)+v(173)*v(204)-v(173)*v(57));
v(166)=-2e0*(-(v(173)*v(190))+v(205)-v(173)*v(58));
v(169)=-2e0*(v(180)*v(204)+v(206)-v(180)*v(57));
S(1)=v(17)*v(181)+mu*(1e0-v(181)/v(14));
S(2)=v(17)*v(174)+mu*(1e0-v(174)/v(14));
S(3)=v(17)*v(177)+mu*(1e0-v(177)/v(14));
S(4)=v(178)*v(207);
S(5)=v(175)*v(207);
S(6)=v(179)*v(207);
S(7)=v(173)*v(207);
S(8)=v(176)*v(207);
S(9)=v(180)*v(207);
dSdE(1,1)=v(181)*v(185);
dSdE(1,2)=v(127);
dSdE(1,3)=v(128);
dSdE(1,4)=v(129);
dSdE(1,5)=v(130);
dSdE(1,6)=v(131);
dSdE(1,7)=v(132);
dSdE(1,8)=v(133);
dSdE(1,9)=v(134);
dSdE(2,1)=v(127);
dSdE(2,2)=v(174)*v(186);
dSdE(2,3)=v(136);
dSdE(2,4)=v(137);
dSdE(2,5)=v(138);
dSdE(2,6)=v(139);
dSdE(2,7)=v(140);
dSdE(2,8)=v(141);
dSdE(2,9)=v(142);
dSdE(3,1)=v(128);
dSdE(3,2)=v(136);
dSdE(3,3)=v(177)*v(189);
dSdE(3,4)=v(144);
dSdE(3,5)=v(145);
dSdE(3,6)=v(146);
dSdE(3,7)=v(147);
dSdE(3,8)=v(148);
dSdE(3,9)=v(149);
dSdE(4,1)=v(133);
dSdE(4,2)=v(141);
dSdE(4,3)=v(148);
dSdE(4,4)=v(150);
dSdE(4,5)=v(151);
dSdE(4,6)=v(152);
dSdE(4,7)=v(153);
dSdE(4,8)=v(191)*(-(mu*v(50))+v(61));
dSdE(4,9)=v(155);
dSdE(5,1)=v(134);
dSdE(5,2)=v(142);
dSdE(5,3)=v(149);
dSdE(5,4)=v(156);
dSdE(5,5)=v(157);
dSdE(5,6)=v(158);
dSdE(5,7)=v(159);
dSdE(5,8)=v(155);
dSdE(5,9)=v(196)*((-mu+v(183))*v(51)+v(175)*v(54));
dSdE(6,1)=v(132);
dSdE(6,2)=v(140);
dSdE(6,3)=v(147);
dSdE(6,4)=v(161);
dSdE(6,5)=v(162);
dSdE(6,6)=v(163);
dSdE(6,7)=v(201)*(-(mu*v(49))+v(60));
dSdE(6,8)=v(153);
dSdE(6,9)=v(159);
dSdE(7,1)=v(131);
dSdE(7,2)=v(139);
dSdE(7,3)=v(146);
dSdE(7,4)=v(165);
dSdE(7,5)=v(166);
dSdE(7,6)=2e0*v(173)*v(197);
dSdE(7,7)=v(163);
dSdE(7,8)=v(152);
dSdE(7,9)=v(158);
dSdE(8,1)=v(129);
dSdE(8,2)=v(137);
dSdE(8,3)=v(144);
dSdE(8,4)=2e0*v(176)*v(200);
dSdE(8,5)=v(169);
dSdE(8,6)=v(165);
dSdE(8,7)=v(161);
dSdE(8,8)=v(150);
dSdE(8,9)=v(156);
dSdE(9,1)=v(130);
dSdE(9,2)=v(138);
dSdE(9,3)=v(145);
dSdE(9,4)=v(169);
dSdE(9,5)=2e0*v(180)*v(202);
dSdE(9,6)=v(166);
dSdE(9,7)=v(162);
dSdE(9,8)=v(151);
dSdE(9,9)=v(157);


function [x]=SMSKDelta(i,j)
if (i==j) , x=1; else x=0; end;
end
function [x]=SMSDeltaPart(a,i,j,k)
l=round(i/j);
if (mod(i,j) ~= 0 | l>k) , x=0; else x=a(l); end;
end
function [x]=Power(a,b)
x=a^b;
end

end