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
% video 16:
%   diff(IR66_video_onsets)(15:18) = [1951,7120,2635];
%   diff(IR68_video_onsets)(15:18) = [1918,7153,2635];
%   diff(CP24_video_onsets)(15:17) = [9245,2627]; NOTE: back 1 because of missing video 16
%   CP24_video_onsets(50) = 61056
%   mean([1951,1918]) = 1935
new(61056+1935:61056+1935+200) = 362.6;
% video 52:
%   diff(IR66_video_onsets)(51:53) = [1584,5370,2634];
%   diff(IR68_video_onsets)(51:53) = [1601,5369,2635];
%   diff(CP24_video_onsets)(50:51) = [7216,2644]; NOTE: back 1 because of missing video 16
%   CP24_video_onsets(50) = 248992
%   mean([1601,1584]) = 1593
new(248992+1593:248992+1593+200) = 362.6;