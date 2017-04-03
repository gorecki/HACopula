function plotbimargins(data, varargin)
%PLOTBIMARGINS - Plots all bivariate projections of data.
%
% Usage:
% 1) plotbimargins(data) - plots all bivariate projections of (copula)-data and show them in 
% unit squares and computes corresponding Kendall's taus
% 2) plotbimargins(data, classes) - plots the data similarly to the
% previous aproach but distiguish them by the provided  classes using
% different colors (7 colors (classes) at maximum)
%
% Example:
% load fisheriris;
% plotbimargins(meas, species);
%
%
% Copyright 2017 Jan Górecki

COLORS = 'bkrgymc';
[n, d] = size(data);

% assure that height = width for the figure
fig = figure;
cla reset;
screenSize = get( 0, 'Screensize' );
maxSize = min(screenSize(3:4)); % get the smaller size of the screen
SIZE_RATIO = 0.8;
maxSize = round(maxSize * SIZE_RATIO); % make it smaller
lowX = round((screenSize(3) - maxSize)/2);
lowY = round((screenSize(4) - maxSize)/2);
set(fig, 'Position', [lowX, lowY, maxSize, maxSize]);


if size(varargin, 1) >= 1
    Y = varargin{1};
    if ~iscell(Y) %to cell transform
        Yc = cell(size(Y,1),1);
        for i = 1:size(Y,1)
            Yc(i) = {num2str(Y(i))};
        end
        Y = Yc;
    end
    noc = unique(Y);
    for k = 1:size(noc,1)
        indices = (strcmp(Y,noc(k)));
        A = corr(data(indices,:), 'type', 'Kendall');
        for i = 1:d
            for j = 1:d
                subplot(d,d,(i - 1)*d + j);
                if (i < j)
                    plot(data(indices,i), data(indices,j),'.','MarkerSize', 5,'MarkerEdgeColor', COLORS(k));
                    hold on;
                    axis([min(data(:,i)) max(data(:,i))  min(data(:,j)) max(data(:,j))]);
                elseif (i == j)
                    axis off;
                    rectangle('Position',[0,0,1,1])
                    text(0.5, 0.5, sprintf('%s%d%s','$U_{',i,'}$'), 'HorizontalAlignment', 'center', 'Interpreter','latex');
                else
                    axis off;
                    rectangle('Position',[0,0,1,1])
                    text(0.5, k/3 , sprintf('%s%d%d%s%3.3f%s','$\tau^n_{',j,i,'}=',A(i,j),'$'), 'HorizontalAlignment', 'center', 'Interpreter','latex', 'FontSize', 8, 'Color', COLORS(k));
                end
            end
        end
    end
    %legend(noc{:});
else
    A = corr(data, 'type', 'Kendall');
    clf;
    ha = tightsubplot(d,d, 0.003);
    for i = 1:d
        for j = 1:d
            axes(ha((i - 1)*d + j));
            if (i < j)
                plot(data(:,i), data(:,j),'o','MarkerSize', 3,'MarkerEdgeColor','k');
                axis off;
                rectangle('Position',[0,0,1,1])
            elseif (i == j)
                axis off;
                rectangle('Position',[0,0,1,1])
                text(0.5, 0.5, sprintf('%s%d%s','$U_{',i,'}$'), 'HorizontalAlignment', 'center', 'Interpreter','latex','FontSize', 11);
            else
                axis off;
                rectangle('Position',[0,0,1,1])
                %text(0.5, 0.5, sprintf('%s%d%d%s%3.3f%s','$\tau^n_{',j,i,'}=',A(i,j),'$'), 'HorizontalAlignment', 'center', 'Interpreter','latex', 'FontSize', 8);
                text(0.5, 0.7, sprintf('%s%d%d%s','$\tau^n_{',j,i,'}=$'), 'HorizontalAlignment', 'center', 'Interpreter','latex', 'FontSize', 11);
                text(0.5, 0.3, sprintf('%s%3.3f%s','$',A(i,j),'$'), 'HorizontalAlignment', 'center', 'Interpreter','latex', 'FontSize', 11);
            end
        end
    end
end


