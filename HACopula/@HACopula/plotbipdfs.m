function fig = plotbipdfs(obj, maxZLim)
%PLOTBIPDFS - plot the probability density function for all bivariate margins
%
% Input (obligatory):
% obj       - A HACopula object.
% Input (optional):
% maxZLim   - Either a positive real value to which the Z-axis of each bivariate
%             plot is truncated to, i.e., zlim([0 maxZLim]) is applied to each
%             bivariate margin plot. Or 'off' if no truncation should be
%             applied.
%           
%
% Output:
% fig    - a handle to the rendered figure
%
%
% Copyright 2018 Jan Gorecki

narginchk(1,2);

if ~exist('maxZLim', 'var')
    % default settings
    useZLim = true;     % truncate the z-axis
    maxZLim = 4;         % truncate the z-axis to 4
else
    if isnumeric(maxZLim)
        useZLim = true;
    elseif strcmp(maxZLim, 'off')
        useZLim = false;
    else
        error('HACopula:plotbipdfs:BadInputs', 'HACopulafit::plotbipdfs: The parameter ''maxZLim'' should be a value  > 0 or ''off''');
    end
end

N_STEPS = 25; % the size of the mesh for each bivariate plot

% set the size of the figure
SIZE_RATIO = 0.8;
fig = figure;
cla reset;
screenSize = get(0, 'Screensize' );
maxSize = min(screenSize(3:4)); % get the smaller size of the screen
maxSize = round(maxSize * SIZE_RATIO); % make it smaller
lowX = round((screenSize(3) - maxSize)/2);
lowY = round((screenSize(4) - maxSize)/2);
set(fig, 'Position', [lowX, lowY, maxSize * 1.15, maxSize]); % 1.15 makes the contour plots perfect square

% remove the axes
ax = gca;
if isoctave
    set(ax, 'visible', 'off');     
  else % MATLAB
    ax.Visible = 'off';
end

d = obj.Dim;
ha = tightsubplot(d, d, [0.01 0.02], [0.01 0.01], [0.01 0.01]);

for i = 1:d
    % plot the bivariate marginal copula densities
    for j = 1:d
        axes(ha((i - 1)*d + j));
        %subplot(d, d, (i-1)*d+j);
        if i < j
            U = linspace(0.05, 0.95, N_STEPS);
            [U1, U2] = meshgrid(U, U);
            biAC = getbimargin(obj, i, j);
            cPdf = ACpdf(biAC.Family, biAC.Parameter, U1(:), U2(:));
            
            % plot above the diagonal
            if isoctave
                surf(U1, U2, reshape(cPdf, N_STEPS, N_STEPS), ...
                    'EdgeColor', 'none');
            else % MATLAB
                surf(U1, U2, reshape(cPdf, N_STEPS, N_STEPS), ...
                    'FaceAlpha' ,0.8, 'EdgeColor', 'none');
            end
            set(gca,'xticklabel',[]);
            set(gca,'yticklabel',[]);
            
            if useZLim
                zlim([0 maxZLim]);
            end
            
            % plot below the diagonal
            axes(ha((j - 1)*d + i));
            contour(U1, U2, reshape(cPdf, N_STEPS, N_STEPS));
            set(gca,'xticklabel',[]);  % remove ticks
            set(gca,'yticklabel',[]);
            hold on;
            
            if isoctave
                text(0.5, 0.5, ['\tau_{' num2str(j) num2str(i) '}='  ...
                    char(10) sprintf('%.3f', biAC.Tau) ], ...
                    'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'FontSize', 13);
            else
                text(0.5, 0.5, ['$\tau_{' num2str(j) num2str(i) '}=$'  ...
                    newline sprintf('%.3f', biAC.Tau) ], ...
                    'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'FontSize', 13);
            end
            
        elseif i == j
            if isoctave
                text(0.5, 0.5, ['U_' num2str(i) ''], 'Interpreter', 'tex', ...
                     'HorizontalAlignment', 'center', 'FontSize', 13);
            else
                text(0.5, 0.5, ['$U_{' num2str(i) '}$'], 'Interpreter', 'latex', ...
                     'HorizontalAlignment', 'center', 'FontSize', 13);
            end
        else
        end
    end
end

end