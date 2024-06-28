function t_d_curve = displacement_curve(amp_displacement,simulation_time,plotflag)
    dt = simulation_time/100;
    dts = dt:dt:simulation_time;
    dts = dts';
    dt_ = dts/simulation_time;
    d = amp_displacement .*  dt_.^3 .* (10 - 15 .* dt_ + 6 .* dt_.^2);
    if plotflag
        figure();
        plot(dts,d,'o');
        title('displacement_curve');
        xlabel('time/s');
        ylabel('d/m')
    end
    t_d_curve = [dts,d];
end