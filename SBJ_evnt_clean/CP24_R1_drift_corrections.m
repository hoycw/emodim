x1 = [256204, 307287];
y1 = [4.3955 87.912];
coeff  = polyfit(x1,y1, 1);
line1  = coeff(1)*[x1(1):x1(2)] + coeff(2);
x2 = [x1(2) 309015];
y2 = [y1(2) -3.5];
coeff2 = polyfit(x2,y2,1);
line2  = coeff2(1)*[x2(1):x2(2)] + coeff2(2);

new = evnt.trial{1};
new(x1(1):x1(2)) = new(x1(1):x1(2))-line1;
new(x2(1):x2(2)) = new(x2(1):x2(2))-line2;

x3 = [313245 315690];
y3 = [-3.4 81.32];
coeff3  = polyfit(x3,y3, 1);
line3  = coeff3(1)*[x3(1):x3(2)] + coeff3(2);
x4 = [x3(2) 316165];
y4 = [y3(2) -2.198];
coeff4 = polyfit(x4,y4,1);
line4  = coeff4(1)*[x4(1):x4(2)] + coeff4(2);

new(x3(1):x3(2)) = new(x3(1):x3(2))-line3;
new(x4(1):x4(2)) = new(x4(1):x4(2))-line4;

% Add in missing videos
% video 16: mean duration in IR66 and IR68 = 19345
%   video 15 onset = 49404
new(49404,49404+19345) = 362.6;
% video 52: mean duration in IR66 anbd IR68 = 