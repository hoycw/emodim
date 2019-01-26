function create_nk_photodiode(subjectm)

% --------------------------------------------------------
% Synchronizes Neuralynx with Nihon Kohden recording, and creates a
% photodiode channel for the latter based on the first
% Tested on IR69 and IR75
% 
% Arjen Stolk, 2018
% astolk@berkeley.edu
% --------------------------------------------------------


% subject-specific settings
chan              = 'RAM2'; % channel used for sync'ing
photo             = 'Photo1'; % Neuralynx photodiode
dset_nl           = '/Users/arjsto/Projects/Ecog/data/IR69/datafiles/2018-02-10_09-24-14'; % Neuralynx dataset
tcgver            = 3; % 1: TCG, 2: TCG kids, 3: dual with IR76 (IR75 only)
d                 = 1; % dataset

% header information
eval(subjectm);

% load nihon kohden
dset_nk           = [subj(tcgver).datadir subj(tcgver).eegfile{d}];
hdr_nk            = ft_read_header(dset_nk);
chanindx_nk       = find(~cellfun(@isempty, regexp(hdr_nk.label, chan))== 1);
coi_nk            = ft_read_data(dset_nk, 'header', hdr_nk, 'chanindx', chanindx_nk);

% load neuralynx
hdr_nl            = ft_read_header(dset_nl);
chanindx_nl       = find(~cellfun(@isempty, regexp(hdr_nl.label, chan))== 1);
coi_nl            = ft_read_data(dset_nl, 'header', hdr_nl, 'chanindx', chanindx_nl);
coi_nl            = coi_nl.*-1; % there might be a flip of sign (should end up with a positive corr)
photoindx_nl      = find(~cellfun(@isempty, regexp(hdr_nl.label, photo))== 1);
if ~isempty(photoindx_nl) % if photodiode is stored in the same folder
  hdr_nl_photo      = hdr_nl;
  photo_nl          = ft_read_data(dset_nl, 'header', hdr_nl, 'chanindx', photoindx_nl);
else % if photodiode is stored separately because of a separate sample frequency (folder_photo)
  hdr_nl_photo      = ft_read_header([dset_nl '_photo']);
  photoindx_nl      = find(~cellfun(@isempty, regexp(hdr_nl_photo.label, photo))== 1);
  photo_nl          = ft_read_data([dset_nl '_photo'], 'header', hdr_nl_photo, 'chanindx', photoindx_nl);
  % additionally correct for timestamp offsets between the two (note,
  % neuralynx thinks at a resolution of 1000kHz, meaning there are 125
  % samples in-between two consecutive samples of a timeseries sampled with
  % 8kHz)
  t_chan = linspace(double(hdr_nl.FirstTimeStamp),double(hdr_nl.LastTimeStamp),hdr_nl.nSamples);
  t_photo = linspace(double(hdr_nl_photo.FirstTimeStamp),double(hdr_nl_photo.LastTimeStamp),hdr_nl_photo.nSamples);
  t_shared = intersect(t_chan,t_photo); 
  if (numel(t_shared)/numel(t_photo))<0.90 % fewer than 90% of samples left
    error('t_shared should be close in numbers to t_photo') % otherwise, there was little overlap in neuralynx sample indices
  end
  i_chan = ismember(t_chan, t_shared);
  i_photo = ismember(t_photo, t_shared);
  coi_nl = coi_nl(find(i_chan==1,1,'first'):find(i_chan==1,1,'last')); % coi_nl(i_chan), but allowing anti-aliasing when downsampling still
  photo_nl = photo_nl(find(i_photo==1,1,'first'):find(i_photo==1,1,'last')); % photo_nl(i_photo)
end
photo_nl          = photo_nl.*-1; % there might be a flip of sign (should end up with a positive corr)

% downsample neuralynx to nihon kohden samplerate
if hdr_nl.Fs>hdr_nk.Fs
  fprintf('downsampling Neuralynx from %d to %d Hz\n', hdr_nl.Fs, hdr_nk.Fs)
  coi_nl_ds         = ft_preproc_resample(coi_nl, hdr_nl.Fs, hdr_nk.Fs, 'resample');
else
  coi_nl_ds         = coi_nl;
end
if hdr_nl_photo.Fs>hdr_nk.Fs
  fprintf('downsampling Neuralynx photodiode from %d to %d Hz\n', hdr_nl_photo.Fs, hdr_nk.Fs)
  photo_nl_ds       = ft_preproc_resample(photo_nl, hdr_nl_photo.Fs, hdr_nk.Fs, 'resample');
else
  photo_nl_ds       = photo_nl;
end

% remove extreme values
coi_nk((coi_nk>median(coi_nk)+4*std(coi_nk))|(coi_nk<median(coi_nk)-4*std(coi_nk))) = median(coi_nk);
coi_nl_ds((coi_nl_ds>median(coi_nl_ds)+4*std(coi_nl_ds))|(coi_nl_ds<median(coi_nl_ds)-4*std(coi_nl_ds))) = median(coi_nl_ds);

% synchronize nihon kohden and neuralynx timeseries
[c, lags]         = xcov(coi_nk', coi_nl_ds');
[~, idx]          = max(c - smooth(c, hdr_nk.Fs*10)); % find the sharp peak

% doublecheck
figure; subplot(2,1,1);
hold on; plot(lags,c);
hold on; plot(lags(idx),c(idx),'k*');
ylabel('correlation');
xlabel('lag');
subplot(2,1,2);
t = 1:numel(coi_nk); 
hold on; plot(t, zscore(coi_nk));
t2 = lags(idx):lags(idx)+numel(coi_nl_ds)-1;
hold on; plot(t2, zscore(coi_nl_ds)+10);
t3 = lags(idx):lags(idx)+numel(photo_nl_ds)-1; % ignore any offset between photo and chan
hold on; plot(t3, zscore(photo_nl_ds)+20);
legend('NK', 'NL', 'NL photo');
print([subj(1).datadir 'datafiles/sync_nk-nl_' subjectm(9:end) '_' num2str(tcgver) '_' num2str(d)], '-dpdf');

% create nihon kohden photodiode channel
trig_ts = ones(1,numel(coi_nk)).*median(photo_nl_ds);
trig_ts(t3(t3>0 & t3<numel(trig_ts))) = photo_nl_ds(t3>0 & t3<numel(trig_ts));
save([subj(tcgver).datadir subj(tcgver).eegfile{d}(1:end-4) '_photo_' num2str(tcgver) '_' num2str(d) '.mat'], 'trig_ts');
