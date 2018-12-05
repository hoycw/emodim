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

%% Add in missing videos
% NEW APPROACH: add in missign videos based on log times + photodiode offset
%   offset is consistent aside from a blip at trial 68
% FORMULA:
%   video_time_n = log_time_n - log_time_1 + photo_time_1 - photo_log_offset
% VIDEO 16:
%   donsets = log_times-video_times; (assuming the videos are matched and intial offsets are removed)
%   mean_offset = mean(donsets([2:15 17:51])); = 0.0624
%   trial_info.log_onset_time(16) = 119.1408;
%   trial_info.log_onset_time(1)  = 67.5962;
%   video_onsets(1) = 11652;
%   CP24_video_onsets(15) = 61056
v16_time = 119.1408 - 67.5962 + 11.652 - 0.0624;
v16_ix = round(v16_time*1000);  % = 63130, which is 2074 past video_onsets(15)
new(v16_ix:v16_ix+200) = 362.6;

% VIDEO 52:
%   trial_info.log_onset_time(52) = 306.8440;
v52_time = 306.8440 - 67.5962 + 11.652 - 0.0624;
v52_ix = round(v52_time*1000);  % = 250834, which is 1842 past video_onsets(51)
new(v52_ix:v52_ix+200) = 362.6;

% ORIGINAL ATTEMPT: Based on matching the event durations from IR66 and IR68
%   Not using the approach because the difference between log and photod
%   onset times in all other CP24 videos is very consistent (~0.0662), and this
%   approach created outliers in that distribution (0.2076 and 0.3168)
% % video 16:
% %   diff(IR66_video_onsets)(15:18) = [1951,7120,2635];
% %   diff(IR68_video_onsets)(15:18) = [1918,7153,2635];
% %   diff(CP24_video_onsets)(15:17) = [9245,2627]; NOTE: back 1 because of missing video 16
% %   CP24_video_onsets(15) = 61056
% %   mean([1951,1918]) = 1935
% new(61056+1935:61056+1935+200) = 362.6;
% % video 52:
% %   diff(IR66_video_onsets)(51:53) = [1584,5370,2634];
% %   diff(IR68_video_onsets)(51:53) = [1601,5369,2635];
% %   diff(CP24_video_onsets)(50:51) = [7216,2644]; NOTE: back 1 because of missing video 16
% %   CP24_video_onsets(50) = 248992
% %   mean([1601,1584]) = 1593
% new(248992+1593:248992+1593+200) = 362.6;