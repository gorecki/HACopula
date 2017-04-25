function fig = plot(obj, varargin)
%PLOT - Visualize a HAC.
%
% Purpose:
% Plots a tree- or a dendrogram-based representation of
% a HACopula object obj.
%
% NOTE:
% 1) The children nodes of each fork are orderer (left to right)
% using the function setlrordering.
% 2) The figure is positioned in the middle of the screen and
% sized 1/2 of the width and height of the screen.
% 3) The figure is plotted in two phases, where the edges are
% plotted during the first phase and the nodes during the second.
% 4) If needed, the description of the forks can be customized
% through genText, whereas the description of the leaves can be
% customized through varText in the auxiliary function
% plotHACrec. Also, the size of the font and of the nodes can be
% customized there.
%
% Inputs (obligatory):
% obj               - a HACopula object
% Inputs (optional):
% GraphType         - A type of the graph ('tree' (default) or
%                     'dendrogram').
% Positioning       - A type of vertical positioning of the forks.
%                     If set to 'even' (default), the forks are
%                     vertically positioned evenly according
%                     to their Level property.
%                     If set to 'tau, the forks are
%                     vertically positioned according
%                     to their Tau property.
% ForkBackground    - If set to 'white' (default), each fork has a white
%                     bar background. If set to 'none', no
%                     background is drawn for a fork.
% MarkerSize        - a positive integer determining the size of the
%                     circle markers used for the leaves
% FontSize          - a positive integer determining the font
%                     size
%
%
% Copyright 2017 Jan Górecki

% default setting of the optional inputs
graphType = 'tree';
positioning = 'even';
forkBackground = 'white';
markerSize = 35; % use 32 for d = 20, 35 otherwise
fontSize = 17; % use 12 for d = 20, 20 otherwise

narginchk(1,7);

% check for additional parameters
if mod(size(varargin,2),2) == 1
    error('HACopula.plot: there must be an even number of the additional parameters');
end

if size(varargin,2) >= 2
    % there are some additional parameters
    parNames = lower(varargin(1:2:size(varargin,2)-1));
    parValues = varargin(2:2:size(varargin,2));
    
    % check of the parameter names are allowed
    for i = 1:size(parNames,2)
        try
            validatestring(parNames{i},{'GraphType', 'Positioning', 'ForkBackground', 'MarkerSize', 'FontSize'});
        catch
            error(['The input, ''' parNames{i} ''', did not match any of the valid parameter names (GraphType, Positioning, ForkBackground, MarkerSize, FontSize).']);
        end
    end
    
    % check for GraphType parameter
    iGraphType = find(strncmp(parNames, lower('GraphType'), 9));
    
    % is GraphType a parameter ?
    if size(iGraphType,2) > 0 %
        if size(iGraphType,2) > 1
            error('HACopula.plot: GraphType is a repeating parameter.')
        end
        graphType = parValues{iGraphType};
        if ~(strcmp(graphType, 'tree') || strcmp(graphType, 'dendrogram'))
            error('HACopula.plot: the value corresponding to GraphType must be ''tree'' or ''dendrogram''.');
        end
    end
    
    % check for Positioning parameter
    iPositioning = find(strncmp(parNames, lower('Positioning'), 11));
    
    % is Positioning a parameter ?
    if size(iPositioning,2) > 0 %
        if size(iPositioning,2) > 1
            error('HACopula.plot: Positioning is a repeating parameter.')
        end
        positioning = parValues{iPositioning};
        if ~(strcmp(positioning, 'even') || strcmp(positioning, 'tau'))
            error('HACopula.plot: the value corresponding to Positioning must be ''even'' or ''tau''.');
        end
        if strcmp(graphType, 'tree') && strcmp(positioning, 'tau')
            warning(['HACopula.plot: Some of the edges could cross each other for the setting GraphType=''tree'' '...
                'and Positioning=''tau''. Use GraphType=''tree'' and Positioning=''even'' instead or use GraphType=''dendrogram''.']);
        end
    end
    
    % check for ForkBackground parameter
    iForkBackground = find(strncmp(parNames, lower('ForkBackground'), 13));
    
    % is ForkBackground a parameter ?
    if size(iForkBackground,2) > 0 %
        if size(iForkBackground,2) > 1
            error('HACopula.plot: ForkBackground is a repeating parameter.')
        end
        forkBackground = parValues{iForkBackground};
        if ~(strcmp(forkBackground, 'white') || strcmp(forkBackground, 'none'))
            error('HACopula.plot: the value corresponding to ForkBackground must be ''white'' or ''none''.');
        end
    end
    
    % check for MarkerSize parameter
    iMarkerSize = find(strncmp(parNames, lower('MarkerSize'), 10));
    % is MarkerSize a parameter ?
    if size(iMarkerSize,2) > 0 %
        if size(iMarkerSize,2) > 1
            error('HACopula.plot: MarkerSize is a repeating parameter.')
        end
        markerSize = parValues{iMarkerSize};
    end
    
    % check for FontSize parameter
    iFontSize = find(strncmp(parNames, lower('FontSize'), 10));
    % is FontSize a parameter ?
    if size(iFontSize,2) > 0 %
        if size(iFontSize,2) > 1
            error('HACopula.plot: FontSize is a repeating parameter.')
        end
        fontSize = parValues{iFontSize};
    end
    
    
end


% get models parameters
maxLevel = getmaxlevel(obj);
d = getdimension(obj);
k = size(obj.Forks,2);

% do a copy of HACopula in order not to change the original
% object by the function setlrordered
cpObj = copy(obj);
% interchange children according to numbers of descendant forks
setlrordering(cpObj);

% start plotting
fig = figure;
cla reset;
screensize = get( 0, 'Screensize' );
set(fig, 'Position', [round(screensize(3)/4), round(screensize(4)/4),...
    round(screensize(3)/2), round(screensize(4)/2)]);
set(gcf,'Units','normal')
%set(gca,'Position',[0.1 0.08 0.88 0.85])
xMarg = 95/screensize(3);
yMarg = 100/screensize(4);
set(gca, 'Position', [xMarg/2 yMarg 1-2*xMarg 1-2*yMarg]);
%set(gca,'ButtonDownFcn','selectmoveresize');

hold on;

% do the recursive plot
plotHACrec(cpObj, 0, maxLevel, 1, d, k, positioning, graphType, forkBackground, 0, markerSize, fontSize);
plotHACrec(cpObj, 0, maxLevel, 0, d, k, positioning, graphType, forkBackground, 0, markerSize, fontSize);
axis off;

%format figure
switch positioning
    case 'even'
        % add a vertical axis showing the levels in a HAC
        %                     axis([0.98 d (-maxLevel - 1.1) -0.7 ]); %[xmin xmax ymin ymax]
        %                     move = 0.1; % shifts horizontally the vertical axis
        %                     set(gca,'XLim',[move - 0.05 d]);
        %
        %                     plot([move move], [-1 -maxLevel], 'k');
        %                     for i = 1:maxLevel
        %                         plot([move (-0.05 + move)], [-i -i], 'k');
        %                         text(move-0.05, -i, num2str(i), 'HorizontalAlignment', 'right', 'FontSize',17, 'Interpreter','latex', 'rotation',0);
        %                     end
        %                     h = text(move-0.5, -(maxLevel-1)/2-1, 'Level', 'HorizontalAlignment', 'center', 'FontSize',17, 'Interpreter','latex', 'rotation',0);
        %                     set(h, 'rotation', 90)
    case 'tau'
        axis([0.98 d 0 1]); %[xmin xmax ymin ymax]
        move = 0.3;
        %set(gca,'YLim',[0 1],'Layer','top');
        set(gca,'XLim',[move - 0.05 d]);
        
        plot([move move], [0 1], 'k');
        plot([move (-0.05 + move)], [0 0], 'k');
        plot([move (-0.05 + move)], [1 1], 'k');
        text(move, 0.5, '$\tau$ ', 'HorizontalAlignment', 'right', 'FontSize',17, 'Interpreter','latex', 'rotation',0);
        text(move, 0, '1 ', 'HorizontalAlignment', 'right', 'FontSize',17, 'Interpreter','latex', 'rotation',0);
        text(move, 1, '0 ', 'HorizontalAlignment', 'right', 'FontSize',17, 'Interpreter','latex', 'rotation',0);
    otherwise
        error(['HACopula::plotHACrec: Positioning ' positioning ' is not supported.']);
end

hold off;

% save the plot to a specified location
%set(fig_handle1, 'PaperPositionMode','auto');
%print(fig_handle1, '-depsc', ['k:\copulas\figures\copula_' ...
%                               num2str(d) 'D_' [ord(:).type] '.eps']);
end