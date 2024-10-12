function motor_data_plot( V, K_v, RPM, RPMLim, RPMIntersec, T, Q, Q_m, I, v, I_max, eta_m, eta)
%MOTOR_DATA_PLOT An interactive plot for different propulsive parameters
	static = false;
	if ~exist('eta','var') %If efficiency data is missing, assume static thrust
		static = true;
	end

	NUM_SUBPLOTS_Y = 4;
	NUM_SUBPLOTS_X = 1;
	
	subplot(NUM_SUBPLOTS_Y,NUM_SUBPLOTS_X,1);
	title(['V = ' num2str(V) ', K_v = ' num2str(K_v)])
	hold on;
	if(~static)
		eta_plot = plot(RPM, eta);
	end
	eta_m_plot = plot(RPM, eta_m);
	hold off;
	xlim(RPMLim);
	ylim([0,1]);
	line([RPMIntersec, RPMIntersec], ylim,'Color','black','LineStyle','--')
	%xlabel('RPM');
	if(~static)
		ylabel('\eta');
		legend([eta_plot,eta_m_plot], {'\eta_{prop}', '\eta_{motor}'},'location','best');
	else
		ylabel('\eta_m');
	end
	
	subplot(NUM_SUBPLOTS_Y,NUM_SUBPLOTS_X,2);
	plot(RPM, T/9.80665);
	xlim(RPMLim);
	%line(xlim,[T_des_kgf,T_des_kgf],'Color','black','LineStyle','--');
	line([RPMIntersec, RPMIntersec], ylim,'Color','black','LineStyle','--')
	%xlabel('RPM');
	ylabel('T (kgf)');

	subplot(NUM_SUBPLOTS_Y,NUM_SUBPLOTS_X,3);
	hold on;
	prop_Q_plot = plot(RPM, Q);
	xlim(RPMLim);
	%xlabel('RPM');
	ylabel('Q (Nm)');

	subplot(NUM_SUBPLOTS_Y,NUM_SUBPLOTS_X,3);
	motor_Q_plot = plot(RPM, Q_m);
	xlim(RPMLim)
	Q_plot_y_lim = ylim;
	ylim([min(0,Q_plot_y_lim(2)-1),Q_plot_y_lim(2)])
	line([RPMIntersec, RPMIntersec], ylim,'Color','black','LineStyle','--')
	legend([prop_Q_plot, motor_Q_plot], {'Q', 'Q_m'},'location','best');
	hold off;

	subplot(NUM_SUBPLOTS_Y,NUM_SUBPLOTS_X,4);
	yyaxis left;
	plot(RPM, I);
	xlim(RPMLim);
	ylim([min(0,1.5*max(I)-1),1.5*max(I)])
	line([RPMIntersec, RPMIntersec], ylim,'Color','black','LineStyle','--')
	line(xlim, [I_max,I_max],'Color','b','LineStyle','--')
	xlabel('RPM');
	ylabel('I (A)');
	yyaxis right;
	plot(RPM,I.*v);
	ylabel('P_{elec} (W)');
	ylim([min(0,1.5*max(I.*v)-1),1.5*max(I.*v)])

end

