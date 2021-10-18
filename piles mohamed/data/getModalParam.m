function [freqs, dampings, other_freqs, lag] = getModalParam(pile, prof, direction, g)
%GETMODALPARAM Summary of this function goes here
%   Detailed explanation goes here
excelFilePath = 'piles mohamed\data\InfosEssaisFrequences.xlsx';

C = readcell(excelFilePath);
M = cell2table(C(2:end, :));
M.Properties.VariableNames = C(1, :);

for kl = 1:size(M, 1) % lines
    if M.pile(kl) == pile && M.profondeur(kl) == prof && M.direction{kl} == direction && M.g(kl) == g
        freqs = {};
        dampings = {};
        kmode = 1;
        while true
            try
                f = M.(sprintf('freq%u', kmode))(kl);
                d = M.(sprintf('damp%u', kmode))(kl);
                
                % freq
                if iscell(f)
                    f = f{1};
                    if ~isnumeric(f)
                        f = str2double(split(f)).';
                    end
                end
                if ismissing(f)
                    f = nan;
                end
                % damping
                if iscell(d)
                    d = d{1};
                    if ~isnumeric(d)
                        d = str2double(split(d)).';
                    end
                end
                if ismissing(d)
                    d = nan;
                end
                
                if any(isnan(f))
                    break
                end
                
                freqs{end+1} = f;
                dampings{end+1} = d;
                kmode = kmode + 1;
            catch
                break
            end
        end
        
        % lag
        lag = M.double_choc_lag(kl);
        if iscell(lag)
            lag = lag{1};
            if ~isnumeric(lag)
                lag = str2double(split(lag)).';
            end
        end
        if ismissing(lag)
            lag = [];
        end
        
        % other freqs
        other_freqs = M.autres_freqs(kl);
        if iscell(other_freqs)
            other_freqs = other_freqs{1};
            if ~isnumeric(other_freqs)
                other_freqs = str2double(split(other_freqs)).';
            end
        end
        if ismissing(other_freqs) % missing or nan
            other_freqs = [];
        end
        
        return
    end
end

error('config not found');

end

