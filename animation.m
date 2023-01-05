function [] = animation(planet,tspan,states,filename)
    %  [] = animation(planet,tspan,states,filename)
    %  
    % saves a .gif file with evolution of states
    % 
    % INPUT:
    %   - planet    :   ASTRO object (to plot the desired planet) i.e. ASTRO("earth")
    %   - times     :   times of the corresponding states [1xN] or [Nx1]
    %   - states    :   [Nx6] matrix containing the states of the orbit [a,e,i,O,w,theta] 
    %   - filename  :   name of the file you wand (no extension!)
    % 
    % --------------------------------------
    % Carlo Zambaldo (info@carlozambaldo.it)
    %       last update: Dec 25, 2022
    % --------------------------------------
    
    scale = 3/2;
    id_fig = figure;
    axis image manual % this ensures that getframe() returns a consistent size
    
    %% pianeta
    %title(filename);
    filename = filename+".gif";
    %axis([-21000 21000 -21000 21000 -21000 21000]./1.35);
    drawnow
    
    % Capture the plot as an image 
    frame_n = 1; % primo frame
    frame = getframe(id_fig); 
    im = frame2im(frame); 
    salva_gif(filename,im,frame_n);
    
    %% orbita iniziale
    for sn = 1:size(states,1)
        %text(-planet.R/2*scale,planet.R/2*scale,planet.R/2*scale,"time: "+num2str(times(sn))./(2*pi/planet.omega)+" [sidereal days]",'HorizontalAlignment','left','FontSize',5000);
        obj = ORBIT.simpleOrbit(states(sn,1),states(sn,2),states(sn,3),states(sn,4),states(sn,5),states(sn,6), planet.mu);
        hold off;
        obj.plotOrbit([0, 2*pi], 0.01, 'r-');
        hold on;
        planet.plotAstro; % disegna la terra
        %view(-130,30);
        xlim(8000.*[-1 1])
        ylim(8000.*[-1 1])
        zlim(8000.*[-1 1])
        title("elapsed time: "+num2str(tspan(sn)./(2*pi/planet.omega))+" [sidereal days]");
        drawnow
        % Capture the plot as an image 
        frame = getframe(id_fig); 
        im = frame2im(frame); 
        salva_gif(filename,im,sn);
    end
end

function [] = salva_gif(filename,im,frame_n)
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if frame_n == 1 % primo frame
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end 
end