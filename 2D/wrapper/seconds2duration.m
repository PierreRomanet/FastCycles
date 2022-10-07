function f = seconds2duration(time_in_seconds)
% Function that convert a duration given in seconds to years, months, days,
% hours, minutes and seconds.

    % Initialisation
    remaining_duration = time_in_seconds;

    
    % Calculate years
    time_in_years = floor(remaining_duration/(365.25*86400));

    % Calculate days
    remaining_duration = time_in_seconds- time_in_years*(365.25*86400);
    time_in_days = floor(remaining_duration/86400);
    
    
    % Calculate hours
    remaining_duration = remaining_duration-86400*time_in_days;
    time_in_hours = floor(remaining_duration/3600);
    
    % Calculate minutes
    remaining_duration = remaining_duration-3600*time_in_hours;
    time_in_minutes = floor(remaining_duration/60);
    
    % Calculate seconds
    remaining_duration = remaining_duration-60*time_in_minutes;
    time_in_seconds =     remaining_duration;
    
    f = [num2str(time_in_years,'%04i') ' years, ' num2str(time_in_days,'%03i') ' days, ' num2str(time_in_hours,'%02i') ' hours, ' num2str(time_in_minutes,'%02i') ' minutes, ' num2str(time_in_seconds,'%06.3f') ' seconds'];
    
end