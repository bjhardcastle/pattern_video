function pattern_video(patpath)%,vidpath,exppath)

% input 1: path to a single pattern file, or a folder of patterns 

% input 2: save path for video output(s) [should be a folder in case of
% multiple video files]
% isValidPatternFile(path)
patPaths = fileorfolder(patpath);

numPats = size(patPaths,1);
for patIdx = 1:numPats
    
    pattern = getPatData(patPaths{patIdx});
    
    
%     getPanelCommands
        % parse .m file
        % default to simple gain from x(1) y(1)
        
   plotvars = getPlotCommands(pattern);
        % flat vs 3d iso
        % blue vs green
        % number of visible panels
        
    plotvars = setupFigure(pattern,plotvars); %adds to existing plotvars
    
    
   for y = 1:pattern.y_num
    vidPath = ['C:\temp\testvid_y=' num2str(y) '.avi'];
    vid = setupVidObj(vidPath);
    open(vid)
    
    for x = 1:pattern.x_num
        updateFrame(x,y,pattern,plotvars)
        
        frame = getframe(plotvars.fig);
        writeVideo(vid,frame)
    end
    
    close(vid)
   end
end

end

function Pats = fileorfolder(path)
% return a cell array of pattern file paths (in full)

% first check whether user-input for pattern file path is a folder or a single file 
    [filepath,name,ext] = fileparts(path);
    
    if isempty(ext)
            % a folder was specified,
            % check the contents contains at least one pattern 
            filelist = what(fullfile(filepath));
            assert(~isempty(filelist.mat),'input folder should contain Pattern_*.mat file(s)')
            
            patfile = zeros(1,length(filelist.mat));
            for pidx = 1:length(filelist.mat)  
                patfile(pidx) = isValidPatternFile(fullfile(filepath,filelist.mat{pidx}));
            end
            assert(any(patfile),'no pattern files found in folder')
            patfilenames = filelist.mat(patfile);
            pathnames = cellstr(repmat(filepath,size(patfilenames,1),1));
            
            Pats = fullfile(pathnames,patfilenames);
            
    elseif ~isempty(ext)
            assert(isValidPatternFile(path),'incorrect file selected')
            % single file was specified
            Pats{1} = path;
    end

end

function tf = isValidPatternFile(path)
    [filepath,name,ext] = fileparts(path);

    tf = 0;
    try
        
    assert(exist(path,'file'),'file doesn''t exist')
        
    % check the naming convention matches a Pattern_*.mat file
    assert(strcmp('.mat',ext),'Pattern file should be a Pattern_*.mat file ')
    assert(strcmp('Pattern_',name(1:8)),'Pattern file should be a Pattern_*.mat file ')
    
    % check the contents of the file (Should contain 'pattern' var)
    matFileData = matfile(path,'Writable',false); % read from .mat file (read-only access to avoid modifying)
    matFileVars = whos(matFileData);
    
    assert(strcmp(matFileVars.name,'pattern'),'Pattern file should containt a ''pattern'' struct')
    %assert(isstruct(matFileVars.pattern),'Pattern file should containt a ''pattern'' struct')

    % everything looks good:
    tf = 1;
    catch ME
        warning(ME.message)
    end
    tf = logical(tf);
end

function pattern = getPatData(patpath)
    % should be called after patpath has been validated as a pattern file
    % (see isValidPatternFile)
    % returns something like:
    %       struct with fields:
    %
    %          gs_val: 4
    %      num_panels: 48
    %           x_num: 2
    %           y_num: 1
    %            Pats: [32×96×2 double]
    %       Panel_map: [4×12 double]
    %     BitMapIndex: [1×48 struct]
    %            data: [3072×1 double]
    
    tempload = load(patpath,'pattern');
    pattern = tempload.pattern;
    
end

function plotvars = getPlotCommands(pattern)
    plotvars = struct;
    plotvars.view = '3d';
    plotvars.panelsColor = 'blue';
    plotvars.FaceAlpha = 0.9;
    plotvars.FaceLighting = 'gouraud';
    plotvars.EdgeColor = 'none';
    plotvars.total_columns = 20; %include panels not visible in arena
    plotvars.visible_panels = ones(size(pattern.Panel_map));
    plotvars.backgroundColor = 0.2*[1 1 1];
    plotvars.axView = [0,25]; % az, el
    
    if strcmp( plotvars.panelsColor, 'blue')
    BlueLEDmap = [ 0.01,0,0.15; ...
                   0.35,0.3,1];
    plotvars.panelsColormap = interp1([0 1],BlueLEDmap,0:1/15:1);
    end
end

function plotvars = setupFigure(pattern,plotvars)
    figHandle = figure('color',plotvars.backgroundColor);
    axHandle = gca;
    
    plotvars.arena_circumference = 8*plotvars.total_columns; % units in 'LED pixels'
    plotvars.arena_radius = 0.5*plotvars.arena_circumference/pi;
    plotvars.arena_height = 8*size(pattern.Panel_map,1);
    
    [x,y,z] = cylinder(plotvars.arena_radius,plotvars.arena_circumference);
    % this returns a one-unit high cylinder.
    % we replicate these arrays along the z- dimension to give an x-unit
    % height cylinder
    for n = 2:plotvars.arena_height+1
        x(n,:) = x(n-1,:);
        y(n,:) = y(n-1,:);
        z(n,:) = z(n-1,:)+1;
    end
    
    % the above arrays map all pixels in the arena, including invisible ones which we don't want to plot
    % modify them so we only plot pattern data to the visible panels
    if plotvars.total_columns ~= size(pattern.Panel_map,2)
        x = x(:,1:8*size(pattern.Panel_map,2)+1);
        y = y(:,1:8*size(pattern.Panel_map,2)+1);
        z = z(:,1:8*size(pattern.Panel_map,2)+1);
        % % [add ability to specify which panels are visible, 
        % % instead of assuming the left-most panels]
    end
    surfHandle = surf(x,y,z);
    daspect([1 1 1])
    pbaspect([1 1 1])
    
    surfHandle.FaceAlpha = plotvars.FaceAlpha;
    surfHandle.EdgeColor = plotvars.EdgeColor;
    surfHandle.FaceLighting = plotvars.FaceLighting;
    figHandle.Colormap = plotvars.panelsColormap;
    
    
    plotvars.surf = surfHandle;
    plotvars.ax = axHandle;
    plotvars.fig = figHandle;
    

    hold( plotvars.ax , 'on') 
    addPanelOutlines(pattern,plotvars)
    addLighting(pattern,plotvars)

end
function addPanelOutlines(pattern,plotvars)
        
    [x2,y2,z2] = cylinder(plotvars.arena_radius,plotvars.total_columns);
    for n = 2:plotvars.arena_height/8+1
        x2(n,:) = x2(n-1,:);
        y2(n,:) = y2(n-1,:);
        z2(n,:) = z2(n-1,:)+8;
    end
    if plotvars.total_columns ~= size(pattern.Panel_map,2)
        x2 = x2(:,1:size(pattern.Panel_map,2)+1);
        y2 = y2(:,1:size(pattern.Panel_map,2)+1);
        z2 = z2(:,1:size(pattern.Panel_map,2)+1);
        % % [add ability to specify which panels are visible, 
        % % instead of assuming the left-most panels]
    end
    
    if ~isfield(plotvars,'surfLines') || isempty(plotvars.surfLines)
        plotvars.surfLines = surf(plotvars.ax,x2,y2,z2);
    end
    
    plotvars.surfLines.FaceColor = 'none';
    plotvars.surfLines.EdgeColor = 0.1.*[1 1 1];

end

function addLighting(pattern,plotvars)

if ~isfield(plotvars,'light') || isempty(plotvars.light)
    plotvars.light = light(plotvars.ax);
end
plotvars.light.Position = [0 0 30];
plotvars.light.Style = 'local';
plotvars.light.Color = plotvars.panelsColormap(end,:);

% l2 = light;
% l2.Position = [0 0 0];
% l2.Color = plotvars.panelsColormap(end,:);
% 
plotvars.surf.FaceLighting = 'gouraud';
plotvars.surf.AmbientStrength = 1;
plotvars.surf.DiffuseStrength = 0.5; 
plotvars.surf.BackFaceLighting = 'reverselit';

plotvars.surf.SpecularStrength = 1;
plotvars.surf.SpecularColorReflectance = 1;
plotvars.surf.SpecularExponent = 7;

end

function updateFrame(x,y,pattern,plotvars)
    frame = flip(pattern.Pats(:,:,x,y),2);
    grayscale = 2^pattern.gs_val;
    framescaled = frame./(grayscale-1);
    
    plotvars.surf.CData = framescaled;
    plotvars.ax.CLim = [0 1];
    plotvars.ax.View = plotvars.axView;
    plotvars.ax.Visible = 'off';
    

end

function vidObj = setupVidObj(vidPath)
    vidObj = VideoWriter(vidPath,'Motion JPEG AVI');
    vidObj.FrameRate = 5;
end