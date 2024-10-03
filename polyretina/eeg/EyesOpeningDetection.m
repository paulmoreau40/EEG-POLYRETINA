EEG = pop_loadset('filename','TestD_EEG_blockT1_1.set','filepath','D:\Data_EEG-AFF\v2\analysis\1_raw-eeglab\TestD');
EEG = pop_resample(EEG, 250);
EEG_filt = custom_filter(EEG, 1, 20);
EEG_sel = pop_select(EEG_filt, 'channel', {'RD1','RR1','R1','Z1','L1','LL1','LD1'});
EEG_epoch = pop_epoch(EEG_sel, {'ObsStart'}, [-10,5]);

figure;topoplot(1:127, EEG.chanlocs, 'style', 'blank', 'electrodes', 'labels');

for tr = 1:EEG_epoch.trials
    pop_plotdata(EEG_epoch, 1, 1:7, tr, sprintf('Trial %d',tr), 1, 1, [0,0]);
end

step = 0.25; % in mV
meanData = squeeze(mean(EEG_epoch.data([2,5,6],:,:),1))*1000; % convert to mV
decMeanData = zeros(size(meanData));
for tr = 1:EEG_epoch.trials
    decMeanData(:,tr) = meanData(:,tr) - step*(tr-1);
end

figure;hold on;
plot(EEG_epoch.times,decMeanData);
xline(0);
for tr = 1:EEG_epoch.trials
    yline(-step*(tr-1));
end
xlabel('Time (ms)');
ylabel('Amplitude (mV)');
