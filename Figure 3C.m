% 1 pre; 2 stim; 3 stim-pre
% 4-D array stim, subj, chans, freqs
addpath(genpath('D:/Studies/REM-REST/Code/EEG/ti_process-main_RR'))
addpath(genpath('D:/Studies/REM-REST/Code/EEG/nonpara_clust_stat'))
local_path = 'D:/Studies/REM-REST/Paper';
diffFreqs = ["lower"]; % 10 Hz difference frequency
clims = {[-2 0]}; % colorbar limits for topoplot
dataTypes = ["log"]; % can be "log" or "z" for spatial z-scored data
load('chanInfo_185') % load channel locations
load("f_ax.mat") % load frequency info
NIter = 1000;
bandNames = ["Sigma"]; % choose power band name
freq_ranges = {[11 16]}; % choose poewr band frequency range
data_name = "NREM";
include = "Y"; % if "Y", exclude protocols non fully in N2
folder = "Stim-Pre_NREM"; % where power data are stored
analyses = ["stim-pre"]; % can be "stim-pre", "stim-post", "post-pre"
save_path = 'C:\Users\sbruno3\OneDrive - UW-Madison\Desktop\Spindle paper';


for i = 1:numel(diffFreqs)
    diffFreq = diffFreqs(i);

    d_path = fullfile(local_path, 'analysis', folder, sprintf('%s', diffFreq));
    d = readtable(fullfile(d_path, sprintf("Freq_%s_metadata.xlsx", diffFreq)), 'Sheet', 'Sheet1');
    if include == "Y"
        included_protos = d.include == "Y";
        d = d(included_protos, :);
    end

    for iDat = 1:numel(dataTypes);
        dataType = dataTypes(iDat);

        if dataType == "log"
            cblab = 'Log Pow Diff';
            cblim = [-1 1];
            cblim_ss = [-2 2];
        end

        if dataType == "z"
            cblab = 'Z-score Pow Diff';
            cblim = [-1 1];
            cblim_ss = [-1 1];
        end

        iFreq = 1;

        for iFreq = 1:numel(freq_ranges)
            psd_path = fullfile(local_path, 'analysis', folder, sprintf('%s', diffFreq));
            data = load(fullfile(psd_path, sprintf('nrem_%s_%s_PSD_xSubjData.mat', dataType, diffFreq)), sprintf('xSubjData_pow_%s', dataType));

            dataa = data.(sprintf('xSubjData_pow_%s', dataType));
            if include == "Y"
                dataa = dataa(:, included_protos, :, :);
            end

            freq_range = freq_ranges{iFreq};
            bandName = bandNames(iFreq);
            freqs = find(f_ax > freq_range(1) & f_ax <= freq_range(2));
            data = squeeze(mean(dataa(:, :, :, freqs), 4));

            groups = findgroups(d(:, 1)); % identify protocols within subjects
            iAn = 1;

            for iAn = 1:numel(analyses)
                analysis = analyses(iAn);
                savestr = sprintf('Cluster_%s_%s_%s_%s', diffFreq, dataType, bandName, analysis);

                % select data of interest form .mat file with EEG power
                if analysis == "stim-pre"
                    d2 = data(1,:,:);
                    d1 = data(2,:,:);
                end

                % remove excess dimension from data
                d1 = squeeze(d1);
                d2 = squeeze(d2);

                % average protocols within subject
                d2 = splitapply(@(x) mean(x, 1), d2, groups);
                d1 = splitapply(@(x) mean(x, 1), d1, groups);

                % initialize params
                nChans = numel(chanInfo);
                alpha = 0.05;

                % obtain N subjects
                if size(d1,1) ~= size(d2,1)
                    error('Number of subjects in d1 and d2 do not match \n')
                end
                NSubj = size(d1,1);

                % calculate the critical value
                % for a two-tailed, paired t-test
                cv = tinv(1-alpha/2, NSubj - 1);

                % for a two-tailed, independent t-test
                % cv = tinv(1-alpha, n1 + n2 - 2);
                fprintf('Performing nonparametric cluster correction topographical stats \n')
                fprintf('Parameters: \n')
                fprintf('alpha: %g, NIter: %g, CV: %g, DoF: %g \n', alpha, NIter, cv, NSubj - 1)

                [contrast, chclus, pval, singleThresholdmax, singleThresholdmin, minclus, pnclus] = cluster_correction_singleThreshold(d1, d2, cv, 'test', 'pt_test', 'nperm', NIter);
                savestr = [savestr '.png'];

                % plot topography PRE, topography STIM, topography of
                % t-statistic of cluster comparison
                figure('units','normalized','outerposition',[0 0 1 1])
                cont_plot = 1;

                subplot(1,3,2)
                topoplot(mean(d1,1), chanInfo, 'maplimits', [clims{i}]);
                cb = colorbar();
                cb.Label.String = 'PSD (log_{10} \muV^2)';
                plot_name = strsplit(analysis, '-');

                title("Sigma Band Power (STIM)")
                fontsize(gcf, 14, 'points')
                set(gca, 'fontweight', 'bold')
                cont_plot = cont_plot + 1;

                subplot(1,3,1)
                topoplot(mean(d2,1), chanInfo, 'maplimits', [clims{i}]);
                cb = colorbar();
                cb.Label.String = 'PSD (log_{10} \muV^2)';
                title("Sigma Band Power (PRE)")
                fontsize(gcf, 14, 'points')
                set(gca, 'fontweight', 'bold')
                cont_plot = cont_plot + 1;

                subplot(1,3,3)
                topoplot(contrast, chanInfo, 'maplimits', [-3 3]);
                % plot significant electrodes, if any
                if size(chclus,2) > 0
                    d3 = zeros(1, nChans);
                    d3(1,chclus) = 0.0005;
                    topoplot(d3, chanInfo, 'style', 'blank', 'plotdisk', 'on', 'headrad', 0.52, 'emarkercolors', {[1 1 1]});
                end

                cb = colorbar();
                title('T-Statistic')
                fontsize(gcf, 14, 'points')
                set(gca, 'fontweight', 'bold')
                cont_plot = cont_plot + 1;

                saveas(gcf, fullfile(save_path, [savestr '.png']))
                saveas(gcf, fullfile(save_path, [savestr '.fig']))
            end
        end
    end
end



