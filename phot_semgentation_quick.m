%% minimal emodim photodiode segmetnation
load CP24_preclean_R1
cfg=[];cfg.channel = {'DC02'};
evnt_ft = ft_selectdata(cfg,data);

[evnt, hdr] = fn_format_data_ft2KLA(evnt_ft);

%% Start SBJ03...
photod_ix = 1;
n_channels = hdr.n_channels;
n_samples = hdr.n_samples;
s_rate = hdr.sample_rate;
data_photo = evnt(photod_ix,:);
data_photo_orig = data_photo;
% data_mic = evnt(mic_ix,:);
% data_mic_orig = data_mic;

% Bring data down to zero
data_photo = data_photo - min(data_photo);

fprintf('\tReading photodiode data\n');
min_event_length = 0.8 * s_rate;    %trial must be at least 0.8 sec (unclear what this is for Alan's videos...)
% Decimate to around 1000Hz (don't worry about aliasing)
if s_rate > 1000
    decimate_v = floor(s_rate/1000);
    data_photo_d = data_photo(1:decimate_v:end);
    min_event_length = floor(min_event_length/decimate_v);
else
    fprintf('================================================\n');
    fprintf('WARNING: evnt sample rate is 1000 or less!!\n');
    data_photo_d = data_photo;
    decimate_v = 1;
end
[~, ~, data_shades] = read_photodiode(data_photo_d, min_event_length, 2);  %2 different shades (video_onset, bsln)
clear data_photo;

data_shades = [diff(data_shades) 0]; % Add a point because diff removes one
video_onsets = find(data_shades>0)'; % 1 to 2,3,4 is word onset. Transpose to make column vector
video_onsets = video_onsets*decimate_v; % convert back to ms
fprintf('\t\tFound %d trials in photodiode channel\n', length(video_onsets));
