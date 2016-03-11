function [H] = couplingLagPlot(R,popts)
%
% Line plot of a coupling statistic (such as I or T) for a single coupling 
% according to lag. 
% 
% Usage: [H] = couplingLagPlot(R,popts)
% 
% ----- Inputs -----
% R = the output structure from a successful run of the ProcessNetwork
%     Software
% popts = a structure of options for plotting. Some are required while
%        others are available as defaults (indicated below). For optional
%        fields, simply don't include the field in the popts structure.
%
% Fields in popts structure (popts.xxx, where xxx is one of the names below)
% NOTE: structure names are case sensitive
% testStatistic = a string indicating a field in the R structure to plot
%       according to lag (this means it must vary according to lag)
% SigThresh (optional) = a string indicating a field in the R structure
%       corresponding with the statistical significance threshold of
%       testStatistic. If absent, no significance threshold will be
%       plotted.
% vars = a 2-element numeric vector or cell array of strings indicating the 
%       variables of the desired coupling to plot. If a cell array, each 
%       string is a variable name within R.varNames. If a numeric vector, 
%       list the two indices of the variables to test within R.varNames. 
%       The order matters here. It will be (1st var > 2nd var), meaning 
%       that the 1st variable corresponds with the row index of the test 
%       statistic and 2nd variable corresponds with the column index.
% fi(optional) = the index of the 4th dimension of the test statistic 
%       corresponding to the fi-th processed file. If absent, the 1st index
%       of the 4th dimension will be used.
% laglim (optional) = a 2-element vector indicating the [min max] lag to
%       show in the plot. If absent, the entire R.lagVect will be shown.
% fignum (optional) = an integer indicating the figure number to plot in.
%       If absent, figure 1 will be used.
% subplot (optional) = a 3-element row vector of integers indicating the
%       [nrows ncolumns index] of the subplot to plot within. This syntax 
%       follows the subplot() function in matlab.
% clearplot (optional) = a binary 0 (no) or 1 (yes) indicating whether to
%       clear the axes before plotting.
% plotProperties (optional) = a cell array with a single row indicating
%       LineSpec properties to give to the plotted test statistic (see 
%       LineSpec in matlab help). The format of this array is such that you
%       simply enclose in curly brackets the inputs that you would normally
%       enter after the x and y variables when using the regular plot 
%       command. For example, if I want to plot the test statistic with 
%       a red dotted line of linewidth 2: 
%       popts.plotProperties = {'linestyle',':','color',[1 0 0],'linewidth',2};
%       If absent, plotting will follow the matlab default
% sigThreshPlotProperties (optional) = same as plotProperties but for the
%       statistical significance line. If absent, the statistical
%       significance threshold will be a dashed line of the same color as
%       the test statistic.
% saveFig (optional) = binary 0 (no) or 1 (yes) to save the figure.
% figName (optional) = a string containing the file name to save the
%       figure, without extension. If absent, then a default name will be
%       used including the name of the test statistic, the file index
%       (popts.fi) and the variable names in the coupling.
% saveFormat (optional) = a string giving the format to save the figure as
%       (without the .). For example: 'png'. See matlab saveas function 
%       for allowable image formats. If absent, the figure will be saved 
%       in the native .fig format.
% outDirectory (optional) = a string giving the full or relative path of
%       the directory to save the figure in. Be sure to include the ending
%       slash (.\myDirectory\). If absent or improper, the figure will be
%       saved in the current directory.
%
% ------ Outputs -------
% The handles of plotted elements are returned in the H structure,
% organized into the following fields:
% fig = the figure handle of the plot
% axis = the axes handle of the plot
% lines = the handle(s) of the plotted lines
% leg = the handle of the legend
%
% --------------------
% Copyright Cove Sturtevant, 2015. 

% Initialize
H.fig = NaN; % figure handle
H.axis = NaN; % axis handle
H.lines = NaN; % line handle
H.leg = NaN; % legend handle

% Check options
if ~isfield(popts,'testStatistic')
    disp('No test statistic indicated. Check options structure.')
    return
elseif ~isfield(R,popts.testStatistic)
    disp([popts.testStatistic ' is not a valid statistic in the results structure.'])
    return
end

% Pull the test statistic
X = eval(['R.' popts.testStatistic]);
[nVars,~,~,nFiles] = size(X);

if ~isfield(popts,'SigThresh')
    XSigThresh = [];
elseif ~isfield(R,popts.SigThresh)
    XSigThresh = [];
else
    XSigThresh = eval(['R.' popts.SigThresh]);    
end

if ~isfield(popts,'vars')
    disp('Missing variable indices of coupling to plot')
    return
elseif numel(popts.vars) < 2
    disp('Two variable indices must be listed in popts.vari')
    return
elseif isnumeric(popts.vars) && (min(popts.vars) < 1 || max(popts.vars) > nVars)
    disp('Invalid variable indices.')
    return
elseif iscell(popts.vars)
    % Get variable indices 
    ri = find(strcmpi(popts.vars(1),R.varNames) == 1,1);
    ci = find(strcmpi(popts.vars(2),R.varNames) == 1,1);
    
    if isempty(ri) || isempty(ci)
        disp('Invalid variable names.')
        return
    end
        
else
    ri = popts.vars(1);
    ci = popts.vars(2);
end
    

if ~isfield(popts,'fi')
    popts.fi = 1;
elseif popts.fi < 1 || popts.fi > nFiles
    disp('Invalid file index (4th dimension).')
    return
end    

if ~isfield(popts,'laglim')
    popts.laglim = [min(R.lagVect) max(R.lagVect)];
elseif isempty(popts.laglim )
    popts.laglim = [min(R.lagVect) max(R.lagVect)];    
elseif numel(popts.laglim) ~= 2
    popts.laglim = [min(R.lagVect) max(R.lagVect)];
elseif popts.laglim(2) - popts.laglim(1) <= 0
    popts.laglim = [min(R.lagVect) max(R.lagVect)];
end    

if ~isfield(popts,'fignum')
    popts.fignum = 1;
elseif isempty(popts.fignum)
    popts.fignum = 1;
elseif ~isnumeric(popts.fignum) || popts.fignum < 0 || isnan(popts.fignum)
    popts.fignum = 1;
end

if ~isfield(popts,'subplot')
    popts.subplot = [];
elseif ~isnumeric(popts.subplot) || numel(popts.subplot) ~= 3 || isnan(sum(popts.subplot))
    disp('Bad subplot option. Plotting as single plot.')
    popts.subplot = [];
end

if ~isfield(popts,'clearplot')
    popts.clearplot = 0;
end

if ~isfield(popts,'saveFig')
    popts.saveFig = 0;
end
if popts.saveFig ~= 0
    if ~isfield(popts,'outDirectory')
        popts.outDirectory = '.\';
    elseif ~exist(popts.outDirectory,'dir')
        disp(['Warning: The output directory ' popts.outDirectory ' does not exist. Setting outDirectory to current directory.'],1);
        popts.outDirectory = '.\';
    end
    
    if ~isfield(popts,'saveFormat')
        popts.saveFormat = 'fig';
    end
end

% Set up figure
figure(popts.fignum); H.fig = gcf;

if ~isempty(popts.subplot)
    subplot(popts.subplot(1),popts.subplot(2),popts.subplot(3))
    H.axis = gca;
    if popts.clearplot ~= 1
        hold on
    else
        delete(H.axis)
        H.axis = subplot(popts.subplot(1),popts.subplot(2),popts.subplot(3));
    end
else
    if popts.clearplot ~= 1
        subplot(1,1,1)
        H.axis = gca;
        hold on
    else
        clf
    end
end

% What lags are we plotting?
lagi = find(R.lagVect >= popts.laglim(1) & R.lagVect <= popts.laglim(2));

% Plot it!
Xp = reshape(X(ri,ci,lagi,popts.fi(1)),[],1);
H.lines = plot(R.lagVect(lagi),Xp);
set(gca,'box','on')

% Modify plotting properties
if isfield(popts,'plotProperties')
    propstr = [];
    for i = 1:length(popts.plotProperties)
        if ischar(popts.plotProperties{i})
            propstr = [propstr '''' popts.plotProperties{i} ''''];
        elseif numel(popts.plotProperties{i}) > 1
            propstr = [propstr '[' num2str(popts.plotProperties{i}) ']'];
        else
            propstr = [propstr num2str(popts.plotProperties{i})];
        end

        if i < length(popts.plotProperties)
            propstr = [propstr ','];
        end
    end
    try
        eval(['set(H.lines,' propstr ')'])
    catch
        disp('Invalid plotProperties option. Plotting with defaults.')
    end
end
    
ylabel([popts.testStatistic])
xlabel('Lag')

% Legend
legText = ['(f' (num2str(popts.fi)) ') ' R.varSymbols{ri} ' > ' R.varSymbols{ci} ]; % legend text for this plot
H.leg = legend(gca); % Get handle whatever legend is here
if ~isempty(H.leg)
    legTextC = H.leg.String;
    legTextC{length(legTextC)+1} = legText;
    legText = legTextC;
end
H.leg = legend(legText);

% Plot statistical significance threshold
if ~isempty(XSigThresh)
    if size(XSigThresh,3) ~= size(X,3)
        XsTp = XSigThresh(ri,ci,1,popts.fi(1))*ones(length(lagi),1);
    else
        XsTp = reshape(XSigThresh(ri,ci,lagi,popts.fi(1)),[],1);
    end
    c = get(H.lines,'color');
    hold on
    H.lines(2,1) = plot(R.lagVect(lagi),XsTp,'--','color',c);
    
    if isfield(popts,'sigThreshPlotProperties')
        propstr = [];
        for i = 1:length(popts.sigThreshPlotProperties)
            if ischar(popts.sigThreshPlotProperties{i})
                propstr = [propstr '''' popts.sigThreshPlotProperties{i} ''''];
            elseif numel(popts.plotProperties{i}) > 1
                propstr = [propstr '[' num2str(popts.sigThreshPlotProperties{i}) ']'];
            else
                propstr = [propstr num2str(popts.sigThreshPlotProperties{i})];
            end

            if i < length(popts.sigThreshPlotProperties)
                propstr = [propstr ','];
            end
        end
        try
            eval(['set(H.lines(2),' propstr ')'])
        catch
            disp('Invalid sigThreshPlotProperties option. Plotting with defaults.')
        end
    end    
        
    % Add to legend
    legText = ['(f' (num2str(popts.fi)) ') ' 'SigThresh ' R.varSymbols{ri} ' > ' R.varSymbols{ci}]; % legend text for this plot
    H.leg = legend(gca); % Get handle whatever legend is here
    if ~isempty(H.leg)
        legTextC = H.leg.String;
        legTextC{length(legTextC)+1} = legText;
        legText = legTextC;
    end
    H.leg = legend(legText);
end
    
if popts.saveFig
    
    if isfield(popts,'figName')
        fname = popts.figName;
    else
        fname = [popts.testStatistic '_' num2str(popts.fi) '_' R.varNames{ri} '_' R.varNames{ci}];
    end
    
    try
        saveas(H.fig,[popts.outDirectory fname],popts.saveFormat);
    catch
        disp(['Problem saving figure: ' [popts.outDirectory fname '.' popts.saveFormat] '. Check options.'])
    end
end

