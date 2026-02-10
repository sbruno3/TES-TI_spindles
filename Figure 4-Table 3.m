% 1 pre; 2 stim; 3 stim-pre
% 4-D array stim, subj, chans, freqs
addpath(genpath('D:/Code/EEG/eeglab2024.1'))
addpath(genpath('D:/Studies/REM-REST/Code/EEG/ti_process-main_RR'))
local_path = 'D:\Studies\REM-REST\Paper\analysis';
titles = {'TES^{15kHz}-TI^{10Hz}' 'TES^{15kHz}-TI^{Peak}' 'TES^{15kHz}' 'No Stimulation'};
diffFreqs = ["lower" "peak" "sham" "baseline"];
load('chanInfo_185.mat') % load channel locations
load('f_ax.mat') % load frequencies info
dataType = ["log"];  % can be "log" or "z" for spatial z-scored data
folder = 'Stim-Pre_NREM';  % where power data are stored
single_subj = "N"; %"Y" to plot individual protocol power spectra
statistic = "npar"; % can be "par" for t-test, "npar" for wilcoxon test, "notest" for no test
bands = {'delta', 'theta', 'sigma', 'beta', 'gamma'}; % choose power band name
bands_range = {[0.5 4] [4 8] [11 16] [16 25] [25 40]}; % choose poewr band frequency range
max_freq = 30; % x axis higher limit for power spectra plots
f_ax = f_ax(f_ax <= max_freq);
include = "Y"; % if "Y", exclude protocols non fully in N2
n_proto = 19; % n of protocol to randomly select for "No Stimulation" analysis

for iFreq = 1:numel(diffFreqs)
    diffFreq = diffFreqs(iFreq);
    d_path = fullfile(local_path, folder, sprintf('%s', diffFreq));
    d = readtable(fullfile(d_path, sprintf("Freq_%s_metadata.xlsx", diffFreq)), 'Sheet', 'Sheet1');
    plot_title = titles{iFreq};
    save_path = fullfile(local_path, folder, diffFreq);

    if include == "Y"
       included_protos = d.include == "Y";
       d = d(included_protos, :);
    end

    groups = findgroups(d(:, 1)); % identify protocols within subjects

    % load all subjects data
    data = load(fullfile(d_path, sprintf('nrem_%s_%s_PSD_xSubjData.mat', dataType, diffFreq)), sprintf('xSubjData_pow_%s', dataType));
    data = data.(sprintf('xSubjData_pow_%s', dataType));
    if include == "Y"
       data = data(:, included_protos, :, :);
    end

    if diffFreq == "baseline"
        idx_sel_protos = [];
        participants = unique(d.sub);
        d_sel = d;

        % selecting one random protocol for each participant
        rng(1234)
        for iPart = 1:numel(participants)
            participantProtocols = d_sel(strcmp(d_sel.sub, participants{iPart}), :);
            rand_num = participantProtocols.abs_proto(randi(numel(participantProtocols.abs_proto)));
            idx_sel_protos = [idx_sel_protos find(strcmp(participants{iPart}, d_sel.sub) & d_sel.abs_proto == rand_num)];
        end

        % Randomly select remaining protocols while balancing 'abs_proto'
        % i.e., making sure that the randomly selected protocols are
        % sampled from, and therefore representative of, the entire nap
        remainingProtocols = d_sel(setdiff(1:height(d_sel), idx_sel_protos), :);

        while numel(idx_sel_protos) < n_proto
            rand_num = randi(height(unique(remainingProtocols.abs_proto)));
            if ismember(rand_num, unique(remainingProtocols.abs_proto))
                proto_pool = remainingProtocols(remainingProtocols.abs_proto == rand_num, :);
                rand_proto = randi(height(proto_pool));
                proto_to_include = proto_pool(rand_proto, :);
                new_proto_ind = find(strcmp(proto_to_include.sub, d_sel.sub) & d_sel.abs_proto == rand_num);
                if ~ismember(new_proto_ind, idx_sel_protos)
                    idx_sel_protos = [idx_sel_protos new_proto_ind];
                end
            end
        end

        d = d_sel(idx_sel_protos, :);
        groups = findgroups(d(:, 1));
        data = data(:,idx_sel_protos,:,:);
    end

    [~, firstIdx] = unique(d.sub, 'first');
    filtered_d = d(firstIdx, [1:3 7:8]);

    % saving single subject data
    for_model = table;
    for_model(1:size(filtered_d, 1), 1:size(filtered_d, 2)) = filtered_d;

    load('f_ax.mat')
    f_ax = f_ax(f_ax <= max_freq);
    data = data(:,:,:,1:find(max_freq <= f_ax, 1));

    pre_subs = squeeze(data(1, :, :, :));
    stim_subs = squeeze(data(2, :, :, :));
    stim_pre_subs = squeeze(data(3, :, :, :));

    if single_subj == "Y"
       savepath_ss = fullfile(save_path, "single_subj");
       if ~exist(savepath_ss)
           mkdir(savepath_ss)
       end
       for iProto = 1:size(d, 1)
           figure; clf; set(gcf, 'Position', [24 74 1747 908])
           pre_ss = squeeze(mean(pre_subs(iProto,:,:), 2));
           stim_ss = squeeze(mean(stim_subs(iProto,:,:), 2));

           plot(f_ax, pre_ss, 'linewidth', 3, 'color', 'k')
           if dataType == "abs"
               set(gca, 'YScale', 'log')
           end
           xlim([0, max_freq])
           hold on
           plot(f_ax, stim_ss, 'linewidth', 3, 'color', 'r')
           legend('pre', 'stim', 'Location', 'southwest')
           savestr = sprintf('Stim-pre_%s_Sub%s_Ses%s_Proto%g', dataType, d.sub{iProto}, d.sess{iProto}, d.abs_proto(iProto));
           title(strrep(savestr, '_', ' '))
           set(gca, 'fontweight', 'bold')
           saveas(gcf, fullfile(savepath_ss, [savestr '.png']))
           saveas(gcf, fullfile(savepath_ss, [savestr '.fig']))

           if size(data, 1) == 6
               figure; clf; set(gcf, 'Position', [24 74 1747 908])
               pre_ss = squeeze(mean(pre_subs(iProto,:,:), 2));
               post_ss = squeeze(mean(post_subs(iProto,:,:), 2));
               plot(f_ax, pre_ss, 'linewidth', 3, 'color', 'k')
               xlim([0, 29.5])
               ylim(ylims_ss)
               hold on
               plot(f_ax, post_ss, 'linewidth', 3, 'color', 'r')
               legend('pre', 'post', 'Location', 'southwest')
               savestr = sprintf('Post-pre_%s_Sub%s_Ses%s_Proto%g', dataType, d.sub{iProto}, d.sess{iProto}, d.abs_proto(iProto));
               title(strrep(savestr, '_', ' '))
               set(gca, 'fontweight', 'bold')
               saveas(gcf, fullfile(savepath_ss, [savestr '.png']))
               saveas(gcf, fullfile(savepath_ss, [savestr '.fig']))


               figure; clf; set(gcf, 'Position', [24 74 1747 908])
               stim_ss = squeeze(mean(stim_subs(iProto,:,:), 2));
               post_ss = squeeze(mean(post_subs(iProto,:,:), 2));
               plot(f_ax, stim_ss, 'linewidth', 3, 'color', 'k')
               xlim([0, 29.5])
               ylim(ylims_ss)
               hold on
               plot(f_ax, post_ss, 'linewidth', 3, 'color', 'r')
               legend('stim', 'post', 'Location', 'southwest')
               savestr = sprintf('Stim-post_%s_Sub%s_Ses%s_Proto%g', dataType, d.sub{iProto}, d.sess{iProto}, d.abs_proto(iProto));
               title(strrep(savestr, '_', ' '))
               set(gca, 'fontweight', 'bold')
               saveas(gcf, fullfile(savepath_ss, [savestr '.png']))
               saveas(gcf, fullfile(savepath_ss, [savestr '.fig']))
           end
       end
    end

    % prepare for grand average
    pre_subs = splitapply(@(x) mean(x, 1), pre_subs, groups);
    pre = squeeze(mean(pre_subs, 1));
        
    stim_subs = splitapply(@(x) mean(x, 1), stim_subs, groups);
    stim = squeeze(mean(stim_subs, 1));

    stim_pre_subs = squeeze(data(3, :, :, :));
    stim_pre_subs = splitapply(@(x) mean(x, 1), stim_pre_subs, groups);
    stim_pre = squeeze(mean(stim_pre_subs, 1));

    p_values = zeros(1, numel(bands));
    median_diff = zeros(1, numel(bands));
    iqr_diff = zeros(1, numel(bands));

    iBand = 1;
    
    % compute median and IQR of the power for each band during STIM and PRE
    % p values of paired test
    for iBand = 1:numel(bands)
        band_range = bands_range{iBand};
        pre_subs_band = pre_subs(:, :, f_ax >= band_range(1) & f_ax < band_range(2));
        pre_band = squeeze(mean(mean(pre_subs_band, 3), 2));
        stim_subs_band = stim_subs(:, :, f_ax >= band_range(1) & f_ax < band_range(2));
        stim_band = squeeze(mean(mean(stim_subs_band, 3), 2));
        if statistic == "par"
           [~, p, ci, stats] = ttest(stim_band - pre_band);
           stat = stats.tstat; 
        end
        if statistic == "npar"
           [p, ~, stats] =  signrank(stim_band - pre_band);
           stat = stats.signedrank;
        end
        p_values(iBand) = p;
        median_diff(iBand) = median(stim_band - pre_band);
        iqr_diff(iBand) = iqr(stim_band - pre_band);

        for_model(:, size(for_model, 2) + 1) = table(pre_band);
        for_model(:, size(for_model, 2) + 1) = table(stim_band);
    end
    
    % save avg stats
    tb = table(bands', median_diff', iqr_diff', p_values')
    filename = sprintf("Table_%s.xlsx", diffFreq);
    writetable(tb,fullfile(save_path, filename));

    % save single subj values
    for_model.Properties.VariableNames = {'sub', 'sess', 'sub_peak', 'freq_stim', 'freq_cat', 'delta_pre', 'delta_stim', 'theta_pre', 'theta_stim', 'sigma_pre', 'sigma_stim', 'beta_pre', 'beta_stim', 'gamma_pre', 'gamma_stim'};
    filename = sprintf("Subj_pow_%s.xlsx", diffFreq);
    writetable(for_model,fullfile(save_path, filename));

    savestr = sprintf('1_stim_pre_%s_PSD_%s_%s', diffFreq, dataType, statistic);

    % plot PRE and STIM power spectra overlapped
    mu_pre = mean(pre);
    mu_stim = mean(stim);

    figure
    plot(f_ax, mu_pre, 'linewidth', 1, 'color', 'k', 'LineStyle',':')
    xlim([0.2 f_ax(end - 2)])

    ylim([-17 20])
    hold on
    plot(f_ax, mu_stim, 'linewidth', 1, 'color', [0.5 0.5 0.5], 'LineStyle','-')
    xlabel("Frequency (Hz)")
    ylabel("PSD (db/Hz)")
    if diffFreq ~= "baseline"
       legendObj = legend('PRE', 'STIM', 'Location', 'northeast');
    end
    if diffFreq == "baseline"
       legendObj = legend('3MIN N2', '3MIN N2', 'Location', 'northeast');
    end

    legendObj.FontSize = 14; % Adjust font size
    legendObj.Position = [0.71, 0.73, 0.15, 0.15]; % [x, y, width, height]
    fontsize(gcf, 14, 'points')
    set(gca, 'fontweight', 'bold')

    saveas(gcf, fullfile(save_path, [savestr '.png']))
    exportgraphics(gcf,fullfile(save_path, [savestr 'band_sig.png']),"Resolution",800)
    saveas(gcf, fullfile(save_path, [savestr '.fig']))
    saveas(gcf, fullfile(save_path, [savestr '.epsc']))

    % plot STIM - PRE difference in power
    % add significance for sigma band, if present
    mu = mean(squeeze(mean(stim_subs, 2)) - squeeze(mean(pre_subs, 2)));
    sd = std(squeeze(mean(stim_subs, 2)) - squeeze(mean(pre_subs, 2)));
    mu = mu(1:numel(mu) - 2);
    sd = sd(1:numel(sd) - 2);
    sd = sd/sqrt(size(pre_subs, 1));
    f_ax = f_ax(1:numel(f_ax) - 2);
    savestr = [plot_title ' PSD diff'];

    figure
    plot(f_ax, mu, 'Color', 'k', 'linewidth', 1.5)
    hold on
    plot(f_ax, mu + sd, 'linewidth', 0.25, 'color', 'k')
    plot(f_ax, mu - sd, 'linewidth', 0.25, 'color', 'k')
    yline(0, 'color', 'r', 'linewidth', 0.5, 'linestyle', '--')
    xlim([0.2 f_ax(end - 2)])
    fontsize(gcf, 16, 'points')
    set(gca, 'fontweight', 'bold')
    xlabel("Frequency(Hz)")
    ylabel("\Delta PSD (db/Hz)")
    fill([f_ax, fliplr(f_ax)], [mu - sd, fliplr(mu + sd)], [0.5 0.5 0.5], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
     if statistic ~= "notest"
       if p_values(3) < 0.05
          plot([11 16], [2 2], 'LineWidth', 1, 'Color', 'r');
          plot(13.5, 2.1, "*", 'Color', 'r');
       end
     end
     title(plot_title)
     ylim([-2 2.9])

    exportgraphics(gcf,fullfile(save_path, [savestr 'band_sig.png']),"Resolution",800)
    saveas(gcf, fullfile(save_path, [savestr 'band_sig.fig']))
end

close all
clear all