function [out1, out2] = plotHACrec(HACModel, iVariable, maxLevel, lines, d, k, positioning, graphType, forkBackground, parentLevel, markerSize, fontSize)
% An auxiliary recursive function for the HACopula's method plot.


% customize the look of the plot
STEP_X = 1; % a horizontal distance of two nodes next to each other
STEP_Y = -1; % % a vertical distance of two nodes from levels i and i+1, respectively

if isa(HACModel, 'HACopula')
    nChildren = size(HACModel.Child,2);
else
    nChildren = 0;
end

%plot the edges
if (nChildren > 1)
    x = zeros(1, nChildren);
    for i = 1:nChildren
        [iVariable, x(i)] = plotHACrec(HACModel.Child{i}, iVariable, maxLevel, lines, d, k, positioning, graphType, forkBackground, HACModel.Level, markerSize, fontSize);
    end
    xm = mean([x(1) x(nChildren)]);
    
    if (lines == 1)
        for i = 1:nChildren
            if isnumeric(HACModel.Child{i})
                % it's a leaf
                switch positioning
                    case 'even'
                        switch graphType
                            case 'dendrogram'
                                plot([x(i) x(i) xm], [(maxLevel + 1) * STEP_Y  HACModel.Level ...
                                    * STEP_Y HACModel.Level * STEP_Y], 'k');
                            case 'tree'
                                plot([x(i) xm], [(HACModel.Level + 1) * STEP_Y  ...
                                    HACModel.Level * STEP_Y], 'k');
                            otherwise
                                error(['HACopula.plotHACrec: GraphType ' graphType ' is not supported.']);
                        end
                    case 'tau'
                        switch graphType
                            case 'dendrogram'
                                plot([x(i) x(i) xm], [0  1-HACModel.Tau ...
                                    1-HACModel.Tau], 'k');
                            case 'tree'
                                plot([x(i) xm], [0  1-HACModel.Tau], 'k');
                            otherwise
                                error(['HACopula.plotHACrec: GraphType ' graphType ' is not supported.']);
                        end
                    otherwise
                        error(['HACopula.plotHACrec: positioning ' positioning ' is not supported.']);
                end
            else
                % it is a fork
                switch positioning
                    case 'even'
                        switch graphType
                            case 'dendrogram'
                                plot([x(i) x(i) xm], [HACModel.Child{i}.Level * STEP_Y ...
                                    HACModel.Level * STEP_Y HACModel.Level * STEP_Y], 'k');
                            case 'tree'
                                plot([x(i) xm], [HACModel.Child{i}.Level * STEP_Y ...
                                    HACModel.Level * STEP_Y], 'k');
                            otherwise
                                error(['HACopula.plotHACrec: GraphType ' graphType ' is not supported.']);
                        end
                    case 'tau'
                        switch graphType
                            case 'dendrogram'
                                plot([x(i) x(i) xm], [1-HACModel.Child{i}.Tau  1-HACModel.Tau ...
                                    1-HACModel.Tau], 'k');
                            case 'tree'
                                plot([x(i) xm], [1-HACModel.Child{i}.Tau 1-HACModel.Tau], 'k');
                            otherwise
                                error(['HACopula.plotHACrec: GraphType ' graphType ' is not supported.']);
                        end
                    otherwise
                        error(['HACopula.plotHACrec: positioning ' positioning ' is not supported.']);
                end
            end
        end
    end
end

%% plot the nodes
if isnumeric(HACModel) % is a leaf
    iVariable = iVariable + 1; %variable counter
    x = iVariable * STEP_X;
    switch positioning
        case 'even'
            switch graphType
                case 'dendrogram'
                    y = (maxLevel + 1) * STEP_Y;
                case 'tree'
                    y = (parentLevel + 1) * STEP_Y;
                otherwise
                    error(['HACopula.plotHACrec: GraphType ' graphType ' is not supported.']);
            end
        case 'tau'
            y = 0;
        otherwise
            error(['HACopula.plotHACrec: positioning ' positioning ' is not supported.']);
    end
    
    if (lines == 0)
        pl = plot(x, y, 'o');
        if isoctave
            set(pl,'MarkerSize',markerSize, 'MarkerEdgeColor','k');
        else % MATLAB
            set(pl,'MarkerSize',markerSize, 'MarkerEdgeColor','k', 'MarkerFaceColor','w');
        end
    end
    %plot_circle(x, y, marker_shift);
else
    x = xm;
    switch positioning
        case 'even'
            y = HACModel.Level * STEP_Y;
        case 'tau'
            y = 1-HACModel.Tau;
        otherwise
            error(['HACopula.plotHACrec: positioning ' positioning ' is not supported.']);
    end
end

if (lines == 0)
    %plot text
    if isnumeric(HACModel)
        type = 'u';
        add = '';
    else
        type = HACModel.Family;
        if type == '?'
            add = '?';
            whiteSpUnknown = '~~~';
        else
            add =  sprintf('%3.3f',HACModel.Parameter);
            whiteSpUnknown = '';
        end
        %add =  ['$\theta_{' num2str(HACModel.TauOrdering) '}$'];
    end
    
    if isnumeric(HACModel)
        if isoctave
            varText = [type '_{', num2str(HACModel) '}'];
            text(x, y, varText, 'HorizontalAlignment', 'center', 'FontSize', fontSize, 'Interpreter', 'tex', 'FontUnits', 'pixels');
        else % MATLAB
            varText = ['$' type '_{', num2str(HACModel), '}$'];
            text(x, y, varText, 'HorizontalAlignment', 'center', 'FontSize', fontSize, 'Interpreter', 'latex', 'FontUnits', 'pixels');
        end
    else
        
        % do centering :)
        if HACModel.TauOrdering < 10
            whiteSpace = '~~~~';
        else
            whiteSpace = '~~~';
        end
        % generators (forks) representation
        % for the outputs for [Górecki et al., 2016b]
        if strcmp(HACModel.Family,'A') || strcmp(HACModel.Family,'C')
            whiteSpaceAdd = '\,';
        else
            whiteSpaceAdd = '~';
        end
        if isoctave
            
            genText = ['\lambda(' num2str(HACModel.TauOrdering) ')' char(10) '(' type ', ' add ')' char(10) '\tau = ' sprintf('%3.3f',HACModel.Tau)];  % [Górecki et al., 2016b] example
            % plot text
            text(x, y, genText,...
                'HorizontalAlignment', 'center', 'FontSize',fontSize, 'Interpreter', 'tex', 'FontUnits', 'pixels',...
                'BackgroundColor', forkBackground, 'Margin', 1);
        else % MATLAB
            genText = ['$' whiteSpace '\lambda(' num2str(HACModel.TauOrdering) ')$' char(10) whiteSpUnknown '(' type ', ' add ')' char(10) whiteSpaceAdd '$\tau = ' sprintf('%3.3f',HACModel.Tau) '$'];  % [Górecki et al., 2016b] example
            % plot text
            text(x, y, genText,...
                'HorizontalAlignment', 'center', 'FontSize',fontSize, 'Interpreter', 'latex', 'FontUnits', 'pixels',...
                'BackgroundColor', forkBackground, 'Margin', 1);
        end
    end
end

out1 = iVariable;
out2 = x;
end