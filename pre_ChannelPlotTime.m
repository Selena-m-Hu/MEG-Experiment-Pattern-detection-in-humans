function pre_ChannelPlotTime(timelock_data, channels_num)
    % This function plots a topography map showing the data for all the MEG
    % sensors using a certain time interval, whose limit values can be set
    % below. It also depicts a butterfly diagram from the temporal data.
    %
    % 
    % timelock_data: data with the structure given after using 
    %   FT_TIMELOCKANALYSIS that will be plotted.
    %
    % channels_num: list of integer indexes indicating the channels that we
    %   are choosing to highlight. They are useful to depict if our automatic
    %   selection algorithm is working properly.
    %
    % Visitor: 
    % Antonio Rodriguez Hidalgo 
    % Dept. of Signal Theory and Communications
    % Universidad Carlos III de Madrid
    % arodh91@gmail.com
    %
    % Principal Investigator:
    % Maria Chait 
    % Ear Institute
    % University College London
    % m.chait@ucl.ac.uk
    %
    % Last update: 07/June/2018
    
    low_t = timelock_data.low_t;
    high_t = timelock_data.high_t;
    
    chns_selected = timelock_data.label(channels_num); 
    
    %% Timelock topography plot
    figure('units','normalized','outerposition',[0 0 1 1]);    
    cfg = [];
    cfg.parameter = 'avg';
    cfg.layout='CTF275.lay';
    cfg.xlim=[low_t, high_t]';
    cfg.marker = 'labels';
    cfg.interactive = 'yes';
    cfg.colorbar = 'yes';

    cfg.markerfontsize = 8;
    cfg.highlight='on';
    cfg.highlightchannel=ft_channelselection(chns_selected, timelock_data.label);
    cfg.highlightfontsize=20;
    subplot(2,6,[1:4,7:10])
    ft_topoplotER(cfg, timelock_data); title ('M100 response');
    
    %% Butterfly diagram
    subplot(2,6,[5:6, 11:12])
    t_chunk = timelock_data.time;
    t_chunk = t_chunk(t_chunk < .3);    
    plot(t_chunk,timelock_data.avg(:,t_chunk < .3)*1E15)
    xlabel('Time (s)')
    ylabel('Magnitude (fT)')
    xlim([min(t_chunk), .3])

    
end
