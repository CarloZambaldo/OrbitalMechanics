function [] = plotGroundTrack(lon, lat, options)
    % [] = plotGroundTrack(lon, lat, options)
    %
    % plots an already propagated ground track, please refer to
    % computeGroundTrack to propagate it.
    %
    % INPUT:
    %   - lon             :  longitude
    %   - lat             :  latitude
    %   - options         :  struct containing all possible options:
    %                         + (Color, LineWidth ...) plotting optional field

if length(lon) ~= length(lat)
    warning("longitude and latitude vectors should be of the same size! Ignoring extra entries.");
end

%% PLOT OPTIONS
    if (nargin < 3) || isempty(options) || ~isfield(options,'LineWidth')
        LineWidth = 1.1;
    else
        LineWidth = options.LineWidth;
    end

    if (nargin < 3) || isempty(options) || ~isfield(options,'Color')
        Color = 'g';
    else
        Color = options.Color;
    end

    if (nargin < 3) || isempty(options) || ~isfield(options,'LineStyle')
        LineStyle = '-';
    else
        LineStyle = options.LineStyle;
    end

    if (nargin < 3) || ~isfield(options,'planet')
        planet = 'Earth';
    else
        planet = options.planet;
    end
    
    %% Draws background of ground track plot.
    % sets background of ground track plot
    if strcmpi(planet,'Earth Coastlines') && (~isfield(options,'planet') || ~strcmpi(options.planet,'none'))
        
        % loads Earth topographic data
        load('topo.mat','topo');
        
        % rearranges Earth topopgrahic data so it goes from -180 to 180 
        % deg longitude
        topoplot = [topo(:,181:360),topo(:,1:180)];
        
        % plots Earth map by making a contour plot of topographic data at
        % elevation of 0
        contour(-180:179,-90:89,topoplot,[0,0],'black');
        
    elseif ~strcmpi(planet,'Blank') && ~strcmpi(planet,'none')
        
        % reads in image
        if strcmpi(planet,'Earth Cloudy')
            cdata = imread('earth.png')+imread('clouds.png');
        elseif strcmpi(planet,'Earth Night')
            cdata = imread('earthnight.png');
        elseif strcmpi(planet,'Earth Night Cloudy')
            cdata = imread('earthnight.png')+0.1*...
                imread('clouds.png');
        else
            cdata = imread(strcat(lower(planet),'.png'));
        end
        
        % plots background
        if ~strcmpi(planet,'none')
            image('CData',flipud(cdata),'XData',[-180,180],'YData',[-90,90]);
        end

        % determines grid style based on background
        if strcmpi(planet,'Sun')
            GridLineWidth = 1;
            GridColor = 'k';
        elseif strcmpi(planet,'Moon')
            GridLineWidth = 0.75;
            GridColor = 'y';
        elseif strcmpi(planet,'Mercury')
            GridLineWidth = 0.75;
            GridColor = 'y';
        elseif strcmpi(planet,'Venus')
            GridLineWidth = 0.5;
            GridColor = 'k';
        elseif strcmpi(planet,'Earth')
            GridLineWidth = 0.5;
            GridColor = 'y';
        elseif strcmpi(planet,'Earth Cloudy')
            GridLineWidth = 1.5;
            GridColor = 'y';
        elseif strcmpi(planet,'Earth Night')
            GridLineWidth = 0.5;
            GridColor = [0.5,0.5,0.5];
        elseif strcmpi(planet,'Earth Night Cloudy')
            GridLineWidth = 0.5;
            GridColor = [0.5,0.5,0.5];
        elseif strcmpi(planet,'Mars')
            GridLineWidth = 0.75;
            GridColor = 'y';
        elseif strcmpi(planet,'Jupiter')
            GridLineWidth = 1;
            GridColor = 'k';
        elseif strcmpi(planet,'Saturn')
            GridLineWidth = 0.65;
            GridColor = 'k';
        elseif strcmpi(planet,'Uranus')
            GridLineWidth = 0.5;
            GridColor = 'k';
        elseif strcmpi(planet,'Neptune')
            GridLineWidth = 0.5;
            GridColor = 'w';
        elseif strcmpi(planet,'Pluto')
            GridLineWidth = 0.65;
            GridColor = 'g';
        end
        
        % add initial and final points
        hold on;
        plot(lon(1),lat(1),'or','LineWidth',1.3);
        plot(lon(end),lat(end),'*r','LineWidth',1.3);

        % manually adds grid
        for i = 1:11
            plot([-180+i*30,-180+i*30],[-90,90],'Color',GridColor, 'LineWidth',GridLineWidth,'LineStyle',':');
        end
        for i = 1:5
            plot([-180,180],[-90+i*30,-90+i*30],'color',GridColor, 'LineWidth',GridLineWidth,'LineStyle',':');
        end
    end
    

    %% Plotting ground track / axis formatting.
    % determines indices where ground track crosses figure border (i.e. from 180 to -180 or -180 to 180) to avoid "jumps" in the plot
    jj = 1;
    for i = 1:length(lon)-1
        if ((lon(i) > 170) && (lon(i+1) < -170)) || ((lon(i) < -170) && (lon(i+1) > 170))
            jj = [jj,i,i+1];
        end
    end
    
    % adds last index to "j" in order to plot ground track between last figure border crossing and the last input longitude
    jj = [jj,length(lon)];
    
    % plots groundtrack (starts new plot every time ground track crosses left border or right border)
    hold on;
    %plot(lon(1:jj(1)), lat(1:jj(1)),'Color',Color,'LineStyle', LineStyle,'LineWidth',LineWidth);
    for i = 2:2:(length(jj))+1
        plot(lon(jj(i-1):jj(i)),lat(jj(i-1):jj(i)), 'Color',Color, 'LineStyle',LineStyle, 'LineWidth',LineWidth, 'HandleVisibility','off');
    end

    % add once again initial and final points
    hold on;
    plot(lon(1),lat(1),'or','LineWidth',1.3);
    plot(lon(end),lat(end),'*r','LineWidth',1.3);
    
    % axis formatting
    axis equal
    grid on
    if strcmpi(planet,'Earth Coastlines') || strcmpi(planet,'Blank')
        ax = gca;
        ax.GridColor = [0.35,0.35,0.35];
        ax.GridAlpha = 1;
        ax.GridLineStyle = ':';
    end

    xlim([-180,180]);
    xticks(-180:30:180);
    ylim([-90,90]);
    yticks(-90:30:90);

    xlabel('Longitude \lambda [deg]', 'FontSize',15);
    ylabel('Latitude \phi [deg]', 'FontSize',15);
    legend("Start Position","End Position");

    title("Ground Track");
    %title("a = "+num2str(obj.a)+"[km], e = "+num2str(obj.e)+"[-], i = "+num2str(obj.i*180/pi)+"[deg], \Omega = "+num2str(obj.O*180/pi)+"[deg], \omega = "+num2str(obj.w*180/pi)+"[deg]");
end