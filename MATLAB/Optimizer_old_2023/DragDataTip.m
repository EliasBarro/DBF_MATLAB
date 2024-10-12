function [txt] = DragDataTip(hFig, varargin)
%DRAGDATATIP makes data tips draggable for most 2-D plots.
% DRAGDATATIP(...) activates the draggable data tips feature for all axes in the current MATLAB
% figure, gcf.
% 
% DRAGDATATIP(HFIG, ...) activates the draggable data tips feature for all axes in the user specifed
% MATLAB figure, HFIG.
% 
% DRAGDATATIP(VARARGIN) activates the draggable data tip features for all axes in the current MATLAB
% figure, gcf and uses the user defined option-pairs defined by VARARGIN.
% 
% DRAGDATATIP(HFIG, VARARGIN) activates the draggable data tip features for all axes in the user
% specifed MATLAB figure, HFIG and uses the user defined option-pairs defined by VARARGIN.
% 
% --------------------------------------------VARIABLES---------------------------------------------
% INPUT VARIABLE(S)
% hFig - Handle to user specified figure. If no inputs are used, then the current figure handle is
%        assigned to this argument. Else if this argument is not a figure handle, then it is treated
%        as the first element of VARARGIN.
% varargin - Provides various feature options using paired input arguments.
%   'draggable' - Choose 'on' (default) or 'off' to create draggable data tips.
%    'fontsize' - A numerical value greater than zero (default = 10). [point units]
%      'header' - Choose 'on' or 'off' (default) to use the legend display name for each data tip.
%                 i.e. - For a data tip added to a line with a display name of 'Filtered Data', the
%                 data tip will include 'Filtered Data' above the typical X- and Y-values.
%   'linestyle' - Choose between line styles to adjust the leader line appearance.
%   'precision' - A string or character vector indicating the display value format. Formats are
%                 similar to those used in fprint() and sprint() for floating point and scientific
%                 value display. (e.g. - '%0.6f' indicates a floating point value with six digits to
%                 the right of the decimal)
%      'xlabel' - A string or character vector. Replaces the default data tip axis label, 'X'.
%      'ylabel' - A string or character vector. Replaces the default data tip axis label, 'Y'.
% 
% OUTPUT VARIABLE(S)
% txt - Cell array of strings used to relabel default MATLAB data tips. Not actually used to return
%       a value to the user.
% 
% -------------------------------------------DESCRIPTION--------------------------------------------
% This function provides a solution to the limited positions allowed by standard MATLAB Data Tips on
% most 2-D plots by replacing them with draggable versions when the "Data Tip" mode is toggled off.
% Since standard data tips can only be moved to a position where one of its four corners touches the
% selected data point, this can unavoidably lead to obstructed portions of graphed data. Draggable
% data tips can be moved anywhere within the current axes, thus allowing better visibility of the
% plot. New visual features are used to reference the draggable data tip to its corresponding data
% point and parent object. These include a connecting line between the data tip textbox and the data
% point as well as a border on the data tip that will try to utilize the same color and line weight
% properties as the target parent object.
% 
% The varargin input argument allows the user to optionally customize the data tips labels. For
% example X- and Y-value labels can be replace with something more descriptive such as 'Time' and
% 'Accel'. In addition to axis labels and when a plot legend is used, a header label can be added to
% the data tip. Headers are added above the axis labels and match the target object's display name
% as it appears in the legend. While the data tip header feature requires a legend to exist prior to
% creation of the data tip, the legend can be deleted after creating the data tip and will not
% effect the header label. If a legend does not exist when the data tip is created, the header
% option is ignored.
% 
% Additional variable input options include setting the numerical precision and/or font size for the
% displayed data values, and adjusting the appearance of the leader lines by specifying the line
% style.
% 
% Note that while the labels/header options are applied to all axes within a single figure, it is
% still possible to use different labels and headers throughout multiple axes or even on a single
% axes for a figure. This is done by applying DRAGDATATIP to a figure and adding a number of new
% draggeble data tips, then reapplying DRAGDATATIP with a new set of options and creating additional
% data tips.
% 
% NOTE: Draggable data tips can currently be removed by using one of the methods listed below:
%   1. Right-click on any single data tip and select delete single or all from the context menu.
%   2. Use the figure's "Edit Plot" tool, select the data tip, and press 'delete'.
%   3. Find and delete axes line objects with the "Tag" property set to 'DraggableDataTip'.
% 
% DRAGDATATIP relies heavily on undocumented features and functions in MATLAB.
% 
% ==================================================================================================
% Special thanks to Francois Bouffard and his draggable() function, for which this function was made
% possible. Necessary portions of draggable() have been incorporate directly to DragDataTip() to
% facilitate bug fixes and addition performance features. I still highly recommend downloading
% draggable() for other uses.
% 
% ---------------------------------------------EXAMPLES---------------------------------------------
% % Create a figure, plot 'peaks', and assign DragDataTip to the figure's axes.
% >>
% hFig = figure;
% plot(peaks);
% grid on
% title('Draggable Data Tips')
% xlabel('Data Index Values')
% ylabel('Acceleration (g''s)')
% DragDataTip
%       OR
% DragDataTip(hFig)
% 
% 
% % Include optional x- and y-axis labels for the data tips, turn on the headers, and specify the
% % display precision for data tip values, font size, and leader line style.
% >>
% hFig = figure;
% hLin = plot(peaks);
% grid on
% title('Draggable Data Tips with Customizable Labels')
% xlabel('Data Index Values')
% ylabel('Acceleration (g''s)')
% lgtxt = arrayfun(@(a) ['Line #',num2str(a)],1:length(hLin),'un',0)';
% hLeg = legend(lgtxt,'Orientation','horizontal','NumColumns',10,'Location','best');
% DragDataTip('xlabel','Index','ylabel','Accel','header','on','fontsize',16,'linestyle',':',...
%   'precision','%0.3g')
%       OR
% DragDataTip(hFig,'xlabel','Index','ylabel','Accel','header','on','fontsize',16,'linestyle',':',...
%   'precision','%0.3g')
% 
% --------------------------------------DEVELOPER INFORMATION---------------------------------------
% AUTHORED BY: Allen Beavers
% DATE AUTHORED: 20-Dec-2018
% 
% CURRENT VERSION: 4.6 - (R2022a) Allen Beavers, 11-Aug-2022
% -----------------------------------------REVISION HISTORY-----------------------------------------
% Version 1.0 - (R2018b) Original version.
% 
% Version 1.1 - (R2018b) Allen Beavers, 21-Dec-2018
%   - Apply custom data tip labels to default MATLAB Data Tips before creating draggable versions.
%   - Fixed code to allow default (gca) assignment of 'hAx' input when not provided.
% 
% Version 1.2 - (R2018b) Allen Beavers, 03-Jan-2019
%   - Updated parsing of varargin.
%   - Fixed bug when applying custom data tip labels to default MATLAB Data Tips.
% 
% Version 2.0 - (R2018b) Allen Beavers, 21-Feb-2019
%   - Changed function to utilize a property listener on the DataCursorManager "Enable" property.
% 
% Version 3.0 - (R2018b) Allen Beavers, 26-Feb-2019
%   - Simplified inputs to allow function usage with a single line of code.
%   - Assigns function to the current figure when the 'hFig' input is not used.
%   - Improved data tip deletion methods.
%   - Fixed bug that allowed data tips to be applied to existing data tip leader lines.
% 
% Version 3.1 - (R2018b) Allen Beavers, 01-Mar-2019
%   - Added optional 'on/off' input to allow header and labels to be customized without creating
%     draggable data tips.
% 
% Version 3.2 - (R2019a) Allen Beavers, 03-Jun-2019
%   - Fixed bug to work on multiple figures simultaneously.
%   - Fixed bug to work with 'duration' and 'calendarduration' data classes.
%   - Added non-draggable date/time text display fixes when using 'datetime' data classes.
%     The "draggable.m" function requires updating for 'datetime' data types to become draggable.
% 
% Version 4.0 - (R2019a) Allen Beavers, 12-Aug-2019
%   - Incorporated snippets of draggable() directly into this function to address multiple issues.
%       - Removed dependency of draggable() in conjunction with this function.
%       - Added functionality to handle 'datetime' data classes.
%       - Fixed bug effecting motion/position of data tips in log scale plots.
%       - Added context menu on right-click to support more data tip deletion options.
%   - Fixed bug to allow use with various other 2-D plots. (e.g. - scatter, area, bar, etc.)
% 
% Version 4.1 - (R2019b) Allen Beavers, 30-Sep-2019
%   - Added option to specify the display precision of data tip values.
%   - Removed the figure's KeyPressFcn dependency.
% 
% Version 4.2 - (R2022a) Allen Beavers, 25-Mar-2022
%   - Updated CreateFcn callback to allow functionality to work without having to reapply
%     DragDataTip to the figure when it is saved and reopened.
%   - Added ability to specify the 'FontSize' value.
%   - Added ability to specify the 'LineStyle' for leader lines.
% 
% Version 4.3 - (R2022a) Allen Beavers, 06-May-2022
%   - Found that DragDataTip was being assigned to the CreateFcn callback during first enabling of
%     the data tip tool and not during function initialization. It is now assigned on initialization
%     of DragDataTip function.
%   - Improved behavior during figure resize, axis limit changes, and adding new plot items.
%   - Sets all opening figures to visible. Programmatically saved figures (not visible) will open
%     visible.
%   - Restructured parts and updated notes to improve performance and readability of code.
% 
% Version 4.4 - (R2022a) Allen Beavers, 20-May-2022
%   - Minor bug fixes and removed unused or redundant code.
%   - Suppressed warning message when converting object to structure.
%   - Data tips are deleted when their associated parent objects are deleted.
%   - Link "Visible" properties between data tip and parent object.
% 
% Version 4.5 - (R2022a) Allen Beavers, 27-May-2022
%   - Bug fix for reinitializing function when opening saved *.fig files.
% 
% Version 5.6 - (R2022a) Allen Beavers, 11-Aug-2022
%   - Removed redundant and/or non-functioning code related to set/get appdata usage.
% 
% ------------------------------REQUIRED USER DEFINED FUNCTIONS/FILES-------------------------------
% DragDataTip.m
% 

%% -----------------------------------------BEGIN FUNCTION------------------------------------------
% Initialize default variables
v = struct('draggable','on','fontsize',10,'header','off','linestyle','-','precision','auto',...
    'xlabel','X','ylabel','Y'); % Optional input defaults
event = []; % Default event variable (required for callback usage)
if nargin==0
    hFig = gcf;
end

% Parse varargin inputs to structure 'v'.
if ~isempty(varargin)
    % Adjust event/varargin inputs
    if ischar(hFig) && isfield(v,lower(hFig))
        varargin = [{hFig},varargin];
        hFig = gcf;
    elseif ~ischar(varargin{1}) || ~isfield(v,lower(varargin{1}))
        event = varargin{1};
        varargin = varargin(2:end);
    end
    p = inputParser;
    validString = @(x) ischar(x) || isstring(x);
    validTxt = @(x) ~isempty(regexp(x,'^%\d*[.]?\d*?[gf]','match','once'));
    validFormat = @(x) validString(x) && (strcmp(x,'auto') || validTxt(x));
    addParameter(p,'draggable',v.draggable,@(x) any(strcmpi(x,{'on','off'})));
    addParameter(p,'fontsize',v.fontsize,@(x) isnumeric(x) && x>0)
    addParameter(p,'header',v.header,@(x) any(strcmpi(x,{'on','off'})));
    addParameter(p,'linestyle',v.linestyle,@(x) ismember(x,{'-','--',':','-.'}));
    addParameter(p,'precision',v.precision,validFormat);
    addParameter(p,'xlabel',v.xlabel,validString);
    addParameter(p,'ylabel',v.ylabel,validString);
    try
        parse(p,varargin{:})
    catch e
        disp(e)
    end
    v = p.Results;
    clearvars p
end

% Determine if function is initializing or if a standard data tip is being created
if nargin<=1 || (isa(hFig,'matlab.ui.Figure') && isempty(event) || ~isprop(event,'EventName'))
    ST = dbstack('-completenames');
    if nargin>1 && isa(event,'matlab.graphics.internal.DataTipEvent') % New standard data tip added
        % Format labeling for standard data tips
        txt = label_DataTip(v,event);
    else % DragDataTip function is being initialized or opened from file
        % Saved figure is being opened
        if isa(hFig,'matlab.ui.Figure') && any(arrayfun(@(a) any(strcmp(a.name, ...
                ["loadFigure","localOpenFigure","openfig"])),ST))
            hFig.Visible = "on";
        end
        % Get current figure handle if not a function input argument
        initialize_func(hFig,v,varargin)
    end    
    return
end

% Get data cursor info
hDCM = datacursormode(gcf); % Data Cursor Manager handle
CI = getCursorInfo(hDCM);

% Create draggable data tips
if ~ishandle(event.Source) && strcmp(hDCM.Enable,'off') && ~isempty(CI) && strcmp(v.draggable,'on')
    % Get the data tip position values
    XPos = arrayfun(@(a) a.Target.XData(a.DataIndex),CI(:));
    YPos = arrayfun(@(a) a.Target.YData(a.DataIndex),CI(:));
    
    % Determine if the axes data is datetime and form new data tip textbox strings accordingly.
    if strcmp(v.precision,'auto')
        Func = @(x) num2str(x);
    else
        Func = @(x) sprintf(v.precision,x);
    end
    if isdatetime(CI(1).Target.XData)
        Fmt = 'dd-mmm-yyyy';
        XStr = arrayfun(@(a) [v.xlabel,': \bf\color{blue}',datestr(a,Fmt),'\color{black}\rm'],...
            XPos,'un',0);
    elseif isduration(XPos(1)) || iscalendarduration(XPos(1))
    	XStr = arrayfun(@(a) [v.xlabel,': \bf\color{blue}',Func(a),'\color{black}\rm'],...
            days(XPos),'un',0);
    else
        XStr = arrayfun(@(a) [v.xlabel,': \bf\color{blue}',Func(a),'\color{black}\rm'],...
            XPos,'un',0);
    end
    if isdatetime(CI(1).Target.YData)
        Fmt = 'dd-mmm-yyyy';
        YStr = arrayfun(@(a) [v.ylabel,': \bf\color{blue}',datestr(a,Fmt),'\color{black}\rm'],...
            YPos,'un',0);
    elseif isduration(YPos(1)) || iscalendarduration(YPos(1))
        YStr = arrayfun(@(a) [v.ylabel,': \bf\color{blue}',Func(a),'\color{black}\rm'],...
            days(YPos),'un',0);
    else
        YStr = arrayfun(@(a) [v.ylabel,': \bf\color{blue}',Func(a),'\color{black}\rm'],...
            YPos,'un',0);
    end
    
    % Adjust textbox header string when the header option is 'on'.
    hAx = gca; % Current axes handle
    if strcmpi(v.header,'on') && ~isempty(hAx.Legend)
        Str = [arrayfun(@(a) ['\bf',a.Target.DisplayName,'\rm'],CI(:),'un',0),XStr,YStr];
    else
        Str = [XStr,YStr];
    end
    
    % Get the current state of the axes next plot property and set to 'add'.
    tmp = hAx.NextPlot;
    hAx.NextPlot = 'add';
    
    % Set the draggable data tip position(s), and create new data tip objects.
    hDT = gobjects([size(Str,1),1]); % Empty obj for new Data tip textbox(es)
    for i=1:size(Str,1)
        XPts = [XPos(i);XPos(i)+diff(hAx.XLim)/100]; % [DataTipXPt; TextObjLeftExtent]
        YPts = YPos(i)*[1;1]; % [DataTipYPt; TextObjLowerExtent]
        hDL = plot(XPts,YPts,v.linestyle,'Tag','DataTipLeader','Marker','o','Color','k',...
            'MarkerFaceColor','k','MarkerSize',3,'PickableParts','none'); % Data tip leader line
        hDT(i) = text(XPts(2),YPts(2),Str(i,:)','Tag','DraggableDataTip','Interpreter','tex',...
                'VerticalAlignment','bottom','BackgroundColor','w','FontSize',v.fontsize);

        % Attempt to match line colors
        Tgt = CI(i).Target;
        try
            cstr = 'matlab.graphics.chart.primitive.'; % Start of class string
            if isa(Tgt,[cstr,'Scatter'])
                if strcmp(Tgt.MarkerEdgeColor,'flat')
                    idx = find(Tgt.XData==CI(i).Position(1) & Tgt.YData==CI(i).Position(2));
                    hDT(i).EdgeColor = Tgt.CData(min([idx,size(Tgt.CData,1)]),:);
                else
                    hDT(i).EdgeColor = Tgt.MarkerEdgeColor;
                end
            elseif any(cellfun(@(x) isa(Tgt,[cstr,x]),{'Area','Bar'}))
                hDT(i).EdgeColor = Tgt.FaceColor;
            else
                hDT(i).EdgeColor = Tgt.Color;
            end
        catch
            hDT(i).EdgeColor = 'k';
        end

        % Assign object userdata, callback functions, and data tip object behavior
        set(hDT(i),'UserData',hDL,'DeleteFcn',@deleteDataTip);
        set(Tgt,'UserData',[Tgt.UserData,hDT(i)],'DeleteFcn',@deleteDataTip);
        hasbehavior(hDL,'legend',false) % Prevent object from being listed in legend
        hasbehavior(hDT(i),'legend',false) % Prevent object from being listed in legend
        
        % Links visiblity properties on objects and saves/appends property listeners to axes appdata
        hlink = linkprop([Tgt,hDT(i),hDL],'Visible');
        setappdata(hAx,'hlink',[hlink;getappdata(hAx,'hlink')])
    end
    
    % Assign axes appdata and draggable data tip callbacks
    arrayfun(@(a) setappdata(a,'xrange',get(hAx,'XLim')),hDT)
    arrayfun(@(a) setappdata(a,'yrange',get(hAx,'YLim')),hDT)
    arrayfun(@(a) set(a,'ButtonDownFcn',@click_object),hDT)
    
    % Append list of handles to the axes appdata property
    setappdata(hAx,'hDT',[hDT;getappdata(hAx,'hDT')]) % Draggable data tip handles
    
    % Delete standard data tips and restore the 'nextplot' axes property to its original state.
    hDCM.removeAllDataCursors(hAx)
    hAx.NextPlot = tmp;
end

end

%% This function is executed when DragDataTip is initially called
function initialize_func(hFig, v, vars)
% hFig - Handle to figure
% v - Optional input defaults
% vars - Initial varargin argument

% Sort variable options structure and reapply to the DragDataTip function and figure handle
persistent hDCMListen

% Clear listeners for deleted objects
if ~isempty(hDCMListen)
    iL = arrayfun(@(a) isvalid(a) & isvalid(a.Object{1}),hDCMListen);
    hDCMListen = hDCMListen(iL);
    nL = length(hDCMListen)+1;
end

% Get Data Cursor Manager handle
hDCM = datacursormode(hFig);

% Setup DataCursorMode listener for 'Enable' property
hFcn = @(src,event)DragDataTip(src,event,...
    'draggable',v.draggable,'fontsize',v.fontsize,'header',v.header,'linestyle',v.linestyle,...
    'precision',v.precision,'xlabel',v.xlabel,'ylabel',v.ylabel);
if isempty(hDCMListen)
    hDCMListen = addlistener(hDCM,'Enable','PostSet',hFcn);
else
    hDCMListen(nL) = addlistener(hDCM,'Enable','PostSet',hFcn);
end

% Set callbacks for Figure, Axes, and Data Cursor Manager
hDCM.UpdateFcn = [{@DragDataTip},vars];
set(hFig,'CreateFcn',[{@DragDataTip},vars],'SizeChangedFcn',@update_extents)
arrayfun(@(a) set(hFig.CurrentAxes.(a),'LimitsChangedFcn',@update_extents),["X","Y"]+"Axis")

end

%% This function is executed when standard MATLAB data tips are added to an axes
function txt = label_DataTip(v, event)
% v - Optional input defaults
% event - Action event handle

% Initialize variables
hAx = gca; % Get current axes handle
pos = cell(2,1);

% Determine whether or not object data is datatime and generate new label text strings
idx = [isdatetime(event.Target.XData),isdatetime(event.Target.YData)];
if strcmp(v.precision,'auto')
    Func = @(x) num2str(x);
else
    Func = @(x) sprintf(v.precision,x);
end
if idx(1)
    pos{1} = datestr(num2ruler(event.Position(1),hAx.XAxis));
else
    pos{1} = Func(event.Position(1));
end
if idx(2)
    pos{2} = datestr(num2ruler(event.Position(2),hAx.YAxis));
else
    pos{2} = Func(event.Position(2));
end
labels = [{v.xlabel};{v.ylabel}];
txt = cellfun(@(x,y) [x,': \bf\color{blue}',y,'\rm\color{black}'],labels,pos,'un',0);

% Update text labels for standard data tips
if strcmpi(v.header,'on') && ~isempty(hAx.Legend) || any(idx)
    if strcmpi(v.header,'on')
        txt = [arrayfun(@(a) ['\bf',a.Target.DisplayName,'\rm'],event,'un',0);txt];
    end
end

end

%% This function is executed when the axis limits of the current axes change
function update_extents(~, ~)
% src - Handle to the action source [not used]
% event - Action event handle [not used]

% Get all draggable data tip handles
hDT = getappdata(gca,'hDT');
if isempty(hDT)
    return
end

% Reorders stack order of plot axes objects
hDT = reorder_stack;

% Update object(s) appdata
arrayfun(@(a) setappdata(a,'xrange',get(gca,'XLim')),hDT)
arrayfun(@(a) setappdata(a,'yrange',get(gca,'YLim')),hDT)
arrayfun(@(a) setappdata(a,'initial_position',a.Position),hDT); % [x,y,w,h]
arrayfun(@(a) setappdata(a,'initial_extent',a.Extent),hDT); % [x,y,w,h]
arrayfun(@(a) setappdata(a,'initial_point',get(gca,'CurrentPoint')),hDT); % [x,y,~]
arrayfun(@(a) movefcn([],[],a),hDT)

end

%% This function is executed when the object is selected via a left or right mouse-click
function click_object(hDT, event)
% hDT - Handle to data tip
% event - Action event handle

% Reorders stack order of plot axes objects
reorder_stack;

% Determine which mouse-button was used
if event.Button==1 % [left-click] Initiate drag motion
    setappdata(hDT,'initial_position',hDT.Position); % [x,y,w,h]
    setappdata(hDT,'initial_extent',hDT.Extent); % [x,y,w,h]
    setappdata(hDT,'initial_point',get(gca,'CurrentPoint')); % [x,y,~]
    set(gcf,'WindowButtonUpFcn',@deactivate_movefcn);
    set(gcf,'WindowButtonMotionFcn',{@movefcn,hDT});
elseif event.Button==3 % [right-click] Create a context menu for the target data tip
    hCM = uicontextmenu;
    hDT.UIContextMenu = hCM;
    uimenu(hCM,'Text','Delete Current Data Tip','Accel','D','MenuSelectedFcn',{@CMSelection,hDT});
    uimenu(hCM,'Text','Delete All Data Tips','Accelerator','A','MenuSelectedFcn',@CMSelection);
end

end

%% This function deletes data tip(s) targeted from the UIContextMenu
function CMSelection(hCM, ~, hDT)
% hCM - Handle to the source context menu
% hDT - Handle to data tip "textbox"

% Delete object based on type
if strcmp(hCM.Label,'Delete Current Data Tip')
    delete(hDT)
elseif strcmp(hCM.Label,'Delete All Data Tips')
    delete(getappdata(gca,'hDT'))
end

end

%% This function deletes user data assigned to an object
function deleteDataTip(src, ~)
% src - Handle to the action source
% event - Action event handle [not used]

% Deletes all data tips and/or corresponding leader lines assigned to the object
arrayfun(@(a) delete(a.UserData),src)

end

%% This function deactivates the WindowButtonMotionFcn for the target figure
function deactivate_movefcn(hFig, ~)
% hFig - Handle for the *figure* containing the object
% event - Action event handle [not used]
% hDT - The draggable object handle

% Reset MotionFcn and ButtonUpFcn callbacks
hFig.WindowButtonMotionFcn = '';
hFig.WindowButtonUpFcn = '';
end

%% This function allows the data tip to be moved to a new position with the parent axes
function movefcn(~, ~, hDT)
% src - Handle to the action source [not used]
% event - Action event handle [not used]
% hDT - Handle to data tip

% Retrieve data saved in the figure and current mouse position
dat = getappdata(hDT);
current_point = get(gca,'CurrentPoint'); % current location of the mouse pointer [x,y,~]
% initial location of the mouse pointer [x,y,~]
% extent of the data tip [x,y,w,h]
% object position in the axes [x,y,~]
% axes limits [xmin, xmax] and [ymin, ymax]

% Determine data class type and convert to numeric
DurFcn = @(a) isduration(a) | iscalendarduration(a);
DTFcn = @(a) isdatetime(a);
if DurFcn(dat.xrange)
    dat.xrange = days(dat.xrange);
elseif DTFcn(dat.xrange)
    dat.xrange = [0,days(diff(dat.xrange))];
end
if DurFcn(dat.yrange)
    dat.yrange = days(dat.yrange);
elseif DTFcn(dat.yrange)
    dat.yrange = [0,days(diff(dat.yrange))];
end

% Linearize log-scale plot values
if strcmp(hDT.Parent.XScale,'log')
    dat.initial_point(1,1) = log10(dat.initial_point(1,1));
    tmp = log10(sum(dat.initial_extent([1,3])));
    dat.initial_extent(1) = log10(dat.initial_extent(1));
    dat.initial_extent(3) = tmp-dat.initial_extent(1);
    dat.initial_position(1,1) = log10(dat.initial_position(1,1));
    dat.xrange(1:2) = log10(dat.xrange(1:2));
    current_point(1,1) = log10(current_point(1,1));
end
if strcmp(hDT.Parent.YScale,'log')
    dat.initial_point(1,2) = log10(dat.initial_point(1,2));
    tmp = log10(sum(dat.initial_extent([2,4])));
    dat.initial_extent(2) = log10(dat.initial_extent(2));
    dat.initial_extent(4) = tmp-dat.initial_extent(2);
    dat.initial_position(1,2) = log10(dat.initial_position(1,2));
    dat.yrange(1:2) = log10(dat.yrange(1:2));
    current_point(1,2) = log10(current_point(1,2));
end

% Compute mouse movement (dpt is [dx,dy])
cpt = current_point(1,1:2);
ipt = dat.initial_point(1,1:2);
dpt = cpt-ipt; % Difference between current and initial position points

% Re-compute new position
newpos = min([max([dat.initial_position(1:2)+dpt;[dat.xrange(1),dat.yrange(1)]]); ...
    [dat.xrange(2),dat.yrange(2)]-dat.initial_extent(3:4)]); % [x,y,~]

% Adjust new positions for log-scale plots
if strcmp(hDT.Parent.XScale,'log')
    newpos(1) = 10.^newpos(1);
end
if strcmp(hDT.Parent.YScale,'log')
    newpos(2) = 10.^newpos(2);
end

% Set the new position which actually moves the object
hDT.Position = newpos;

% Get data point value and extents for its corresponding data tip textbox
XPts = hDT.UserData.XData;
YPts = hDT.UserData.YData;
Pos = hDT.Extent;
Pts = Pos(1:2)+[[0,0];Pos(3:4);[Pos(3),0];[0,Pos(4)]];

% Determine if either the x- or y-axis are datetime or duration data types, then calculate the
% distance between the data tip point and each corner of the textbox using the appropriate method.
DurFcn = @(s) isduration(s) | iscalendarduration(s);
DTFcn = @(s) isdatetime(s);
xpts = XPts;
ypts = YPts;
if DurFcn(XPts)
    xpts = days(XPts);
elseif DTFcn(XPts)
    XLim = getappdata(hDT,'xrange');
    xpts = days(XPts-XLim(1));
end
if DurFcn(YPts)
    ypts = days(YPts);
elseif DTFcn(YPts)
    YLim = getappdata(hDT,'yrange');
    ypts = days(YPts-YLim(1));
end
Dis = sqrt((Pts(:,1)-xpts(1)).^2+(Pts(:,2)-ypts(1)).^2);

% Update the end of the leader line value to match the nearest textbox corner position.
if DTFcn(XPts)
    xpos = [XPts(1),Pts(min(Dis)==Dis,1)+XLim(1)];
else
    xpos = [XPts(1),Pts(min(Dis)==Dis,1)];
end
if DTFcn(YPts)
    ypos = [YPts(1),Pts(min(Dis)==Dis,1)+YLim(1)];
else
    ypos = [YPts(1),Pts(min(Dis)==Dis,2)];
end
hDT.UserData.XData = xpos;
hDT.UserData.YData = ypos;

end

%% This function reorders the stacking order of plot items
function hDT = reorder_stack()
% hDT - Handle to all text objects within the current axes

% Initialize variables
hDT = getappdata(gca,'hDT'); % Data tip textbox handle(s)
idx = isvalid(hDT);
if any(~idx)
    % Removes handles for deleted objects from list
    setappdata(gca,'hDT',hDT(idx))
    hDT = getappdata(gca,'hDT');
end
hDL = arrayfun(@(a) a.UserData,hDT); % Data tip leader line handle(s)
hCL = findall(gca,'Type','ConstantLine'); % Constant line handle(s). Requires special treatment
warnID = warning('off','MATLAB:structOnObject');
sCL = arrayfun(@(a)struct(a),hCL); % Get undocumented properties
warning(warnID)

% Rearrange the stack order of the data tip objects.
arrayfun(@(a) set(a.Edge,'Layer','middle'),sCL); % Set ConstantLine uistack
arrayfun(@(a) uistack(a,'top'),hDL)
arrayfun(@(a) uistack(a,'top'),hDT)

end

% -------------------------------------------END FUNCTION-------------------------------------------
% -------------------------------------Written by Allen Beavers-------------------------------------

