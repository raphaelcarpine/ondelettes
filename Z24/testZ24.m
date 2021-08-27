clear all

measure = 1;
config = 1;

[X, t, labels] = getDataZ24(measure, config);

%% test sensor position

defaultNbCh = 33;

disp('channels test:');
for config = 1:9
    [~, ~, labels0] = getDataZ24(1, config);
    flag = true;
    for measure = 2:17
        [~, ~, labels] = getDataZ24(measure, config);
        if ~isequal(labels, labels0)
            flag = false;
            break
        end
    end
    if flag
        fprintf('config #%u: OK', config);
        fprintf((length(labels0) ~= defaultNbCh) + 1, ' (%u channels)\n', length(labels0));
    else
        fprintf('config #%u: ', config);
        fprintf(2, 'problem\n');
        
        groups = {1};
        groupsNbCh = {length(labels0)};
        for measure = 2:17
            [~, ~, labels] = getDataZ24(measure, config);
            flag = false;
            for kg = 1:length(groups)
                [~, ~, labels0] = getDataZ24(groups{kg}(1), config);
                if isequal(labels, labels0)
                    groups{kg}(end+1) = measure;
                    flag = true;
                    break
                end
            end
            if ~flag
                groups{end+1} = measure;
                groupsNbCh{end+1} = length(labels);
            end
        end
        disp('groups:');
        for kg = 1:length(groups)
            fprintf('%u ', groups{kg});
            fprintf((groupsNbCh{kg} ~= defaultNbCh) + 1, ' (%u channels)\n', groupsNbCh{kg});
        end
    end
end


%% test nb points

defaultNbPt = 8*8192;

fprintf('\nnb points test:\n');
for measure = 1:17
    fprintf('measure #%02u: ', measure);
    for config = 1:9
        X = getDataZ24(measure, config);
        fprintf((size(X, 2) ~= defaultNbPt) + 1, '%u ', size(X, 2));
        
        if any(isnan(X), 'all')
            fprintf(2, '(NaNs) ');
        end
        
%         nZeroPad = 25; % zero padding, n points minimum test
%         X0 = true(size(X) - [0, nZeroPad-1]);
%         for kz = 1:nZeroPad
%             X0 = X0 & X(:, kz:end-nZeroPad+kz) == 0;
%         end
%         if any(X0, 'all')
%             fprintf(2, '(0s) ');
%         end
        
        
        fprintf('; ');
    end
    fprintf('\b\b\n');
end


%% plot

while true
    fprintf('\n');
    cont = input('plot signal? ', 's');
    if ~ismember({cont}, {'', 'y', 'yes'})
        break
    end
    
    while true
        try
            measure = input('measure: ');
            config = input('config: ');
            if measure == 0 && config ~= 0
                for measure = 1:17
                    [X, ~, labels] = getDataZ24(measure, config);
                    fig = figure('Name', sprintf('measure %02u, config %u', [measure, config]));
                    plot(X.');
                    legend(labels, 'NumColumns', 2);
                    selectLine(gca);
                end
                break
            elseif config == 0
                for config = 1:9
                    [X, ~, labels] = getDataZ24(measure, config);
                    fig = figure('Name', sprintf('measure %02u, config %u', [measure, config]));
                    plot(X.');
                    legend(labels, 'NumColumns', 2);
                    selectLine(gca);
                end
                break
            else
                [X, ~, labels] = getDataZ24(measure, config);
                fig = figure('Name', sprintf('measure %02u, config %u', [measure, config]));
                plot(X.');
                legend(labels, 'NumColumns', 2);
                selectLine(gca);
                break
            end
        catch
            continue
        end
        break
    end
end