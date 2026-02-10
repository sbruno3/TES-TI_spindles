addpath(genpath('D:/Code/EEG/eeglab2024.1'))
addpath(genpath('D:/Studies/REM-REST/Code/EEG/ti_process-main_RR'))

root_path = 'Z:/REM-REST/spindle_detection/spindles/data';
results_folder_path = 'Z:/REM-REST/spindle_detection/spindles/results';
which_protos_file = 'Z:/REM-REST/spindle_detection/spindles/Spindles_protos.xlsx';
save_path = fullfile(results_folder_path, 'int_spindle_activity');
if ~exist(save_path)
    mkdir(save_path)
end
spindle_range = [11 16];

d = readtable(which_protos_file);
protos_table = d(:, [1:6]);

for iProto = 1:size(d,1)
    % load EEG and filter in the spindle range
    sub = sprintf('%03d', d.Sub(iProto));
    sess = d.Sess{iProto};
    abs_proto = d.abs_proto(iProto);

    EEG = pop_loadset(fullfile(root_path, sub, sess, sprintf('REMREST_%s_%s_forYasa.set', sub, sess)));
    EEG = pop_eegfiltnew(EEG, spindle_range(1), spindle_range(2), [], 0, [], 0);
    
    % find pre and stim intervals and define intervals duration
    % ramping periods excluded
    for iEv = 1:size(EEG.event, 2)
        if abs_proto == EEG.event(iEv).proto_ind
            if strcmp(EEG.event(iEv).type, 'pre start')
                pre_start = EEG.event(iEv).latency;
            end
            if strcmp(EEG.event(iEv).type, 'stim start') % ramp up start
                pre_end = EEG.event(iEv).latency;
            end
            if strcmp(EEG.event(iEv).type, 'max stim start') % ramp up end
                stim_start = EEG.event(iEv).latency;
            end
            if strcmp(EEG.event(iEv).type, 'max stim end') % ramp down start
                stim_end = EEG.event(iEv).latency;
            end
        end
    end

    pre_dur = pre_end - pre_start;
    stim_dur = stim_end - stim_start;


    % load list of detected spindles
    spindles_list = readtable(fullfile(results_folder_path, sprintf('%s_%s_spindles.csv', sub, sess)));

    % for each channel in ascending order
    chans = unique(spindles_list.Channel);

    q = regexp(chans, '\d+', 'match');
    q = cellfun(@str2double, q);
    [~, idx] = sort(q);
    chans = chans(idx);

    int_activity_table = table('Size', [numel(chans) 7], 'VariableTypes', {'string','double','double','double','double','double','double'}, ...
        'VariableNames', {'Chan','Spin_Int_Act_PRE','Spin_Int_Act_STIM', 'N_Spin_PRE', 'N_Spin_STIM', 'Dur_PRE', 'Dur_STIM'});
    iChan = 1;

    for iChan = 1:numel(chans)
        chan = chans{iChan};
        for iChanLabel = 1:numel(EEG.chanlocs)
            if strcmp(EEG.chanlocs(iChanLabel).labels, chan)
                chan_idx = iChanLabel;
            end
        end

        % rectify channel
        data = abs(EEG.data(chan_idx, :));

        % find and integrate spindles in pre and stim
        int_activity_pre = 0; int_activity_stim = 0; spin_pre = 0; spin_stim = 0;
        spindles_chan = spindles_list(strcmp(spindles_list.Channel, chan), :);

        for iSpin = 1:size(spindles_chan, 1)
            if spindles_chan.Start(iSpin)*EEG.srate >= pre_start & spindles_chan.End(iSpin)*EEG.srate <= pre_end
                int_activity_pre = int_activity_pre + sum(data(pre_start:pre_end));
                spin_pre = spin_pre + 1;
            end
            if spindles_chan.Start(iSpin)*EEG.srate >= stim_start & spindles_chan.End(iSpin)*EEG.srate <= stim_end
                int_activity_stim = int_activity_stim + sum(data(stim_start:stim_end));
                spin_stim = spin_stim + 1;
            end
        end
        int_activity_pre = int_activity_pre/pre_dur;
        int_activity_stim = int_activity_stim/stim_dur;

        % create one Excel file for each participant with integrated
        % spindle activity in pre and stim by channel
        int_activity_table.Chan(iChan) = chan;
        int_activity_table.Spin_Int_Act_PRE(iChan) = int_activity_pre;
        int_activity_table.Spin_Int_Act_STIM(iChan) = int_activity_stim;
        int_activity_table.N_Spin_PRE(iChan) = spin_pre;
        int_activity_table.N_Spin_STIM(iChan) = spin_stim;
        int_activity_table.Dur_PRE(iChan) = pre_dur;
        int_activity_table.Dur_STIM(iChan) = stim_dur;
    end

    table_name = sprintf('%s_%s_proto%g_int_spin_act.csv', sub, sess, abs_proto)
    writetable(int_activity_table, fullfile(save_path, table_name))

    protos_table.mean_int_spin_act_PRE(iProto) = mean(int_activity_table.Spin_Int_Act_PRE);
    protos_table.mean_int_spin_act_STIM(iProto) = mean(int_activity_table.Spin_Int_Act_STIM);

    protos_table.sum_int_spin_act_PRE(iProto) = sum(int_activity_table.Spin_Int_Act_PRE);
    protos_table.sum_int_spin_act_STIM(iProto) = sum(int_activity_table.Spin_Int_Act_STIM);
end

big_table_name = 'Int_spin_act_summary.csv';
writetable(protos_table, fullfile(save_path, big_table_name))





