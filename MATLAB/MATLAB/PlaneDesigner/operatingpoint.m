classdef operatingpoint
    %UNTITLED8 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        v_inf=20;
        alpha=0;
        beta=0;
        p=0;
        q=0;
        r=0;
        flap=0;
        aileron=0;
        elevator=0;
        rudder=0;
%         alpha='alpha = 0.0';
%         beta='beta = 0';
%         p='pb/2V = 0';
%         q='qc/2V = 0';
%         r='rb/2V = 0';
%         flap='flap = 0';
%         aileron='Cl roll mom = 0';
%         elevator='Cm pitchmom = 0';
%         rudder='Cn yaw mom = 0';

    end
    
    methods
        function obj = operatingpoint(v_inf,alpha,beta)
            obj.v_inf=v_inf;
            obj.alpha=alpha;
            obj.beta=beta;
        end
%         function obj = writeRUN(obj,filename)
%             filename=sprintf('%s.run',filename);
%             
%             fh=fopen(filename,'w');
%             
%             %fprintf(fh,'\n------------\nRun case 1: -runcase-\n\n');
%             fprintf(fh,'\n ---------------------------------------------\n Run case  1:   MATLAB-written Run \n\n');
%             
%             fprintf(fh,' alpha        ->  alpha          =   1.16647    \n');
%             fprintf(fh,' beta -> %s\n',obj.beta);
%             fprintf(fh,' p -> %s\n',obj.p);
%             fprintf(fh,' q -> %s\n',obj.q);
%             fprintf(fh,' r -> %s\n',obj.r);
%             fprintf(fh,' flap -> %s\n',obj.flap);
%             fprintf(fh,' aileron -> %s\n',obj.aileron);
%             fprintf(fh,' elevator -> %s\n',obj.elevator);
%             fprintf(fh,' rudder -> %s\n',obj.rudder);
%             
%             fclose(fh);
%         end
        function obj = setalpha(obj,value)
            obj.alpha=sprintf('alpha = %f',value);
        end
    end
end

