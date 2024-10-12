classdef avlplane
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name='tempplane';
        mach=0;
        Sref=1;
        Cref=1;
        Bref=1;
        XYZref=[0 0 0];
        CDp=0; %default profile drag coefficient
        
        surfaces=[];
    end
    
    methods
        function obj=avlplane(surfaces,Sref,Cref,Bref,XYZref)
            obj.surfaces=surfaces;
            if nargin>1
                obj.Sref=Sref;
                obj.Cref=Cref;
                obj.Bref=Bref;
                obj.XYZref=XYZref;
            end
        end
            
        function []=writeAVL(obj,filename)
            filename=sprintf('%s.avl',filename);
            
            fh=fopen(filename,'w');
            
            fprintf(fh,'%s | case title\n',obj.name);
            fprintf(fh,'%f | Mach\n',obj.mach);
            fprintf(fh,'0 0 0 | iYsym iZsym Zsym\n');
            fprintf(fh,'%f %f %f | Sref Cref Bref\n',obj.Sref,obj.Cref,obj.Bref);
            fprintf(fh,'%f %f %f | Xref Yref Zref\n',obj.XYZref(1),obj.XYZref(2),obj.XYZref(3));
            fprintf(fh,'%f | CDp (optional)\n',obj.CDp);

            
            for i=1:length(obj.surfaces)
                s=obj.surfaces(i);
                fprintf(fh,'\n#============\n');
                fprintf(fh,'SURFACE\n');
                fprintf(fh,'%s | surface name\n',s.name);
                fprintf(fh,'%i %f %i %f | Nchord Cspace Nspan Sspace\n',s.Nchord,s.Cspace,s.Nspan,s.Sspace);
                fprintf(fh,'YDUPLICATE\n');
                fprintf(fh,'%i\n',s.Yduplicate);
                fprintf(fh,'ANGLE\n');
                fprintf(fh,'%f\n',s.angle);
                
                for j=1:length(s.sections)
                    sect=s.sections(j);
                    fprintf(fh,'#------------\n');
                    fprintf(fh,'SECTION\n');
                    fprintf(fh,'%f %f %f %f %f %f %f | Xle Yle Zle Chord Ainc Nspan Sspace\n',sect.XYZle(1),sect.XYZle(2),sect.XYZle(3),sect.chord,sect.ainc,sect.Nspan,sect.Sspace);
                    
                    if contains(sect.airfoil,'NACA','IgnoreCase',true)
                        NACAnumberstring=extractAfter(sect.airfoil,'NACA');
                        fprintf(fh,'NACA\n');
                        fprintf(fh,'%s\n',NACAnumberstring);
                    else
                        fprintf(fh,'AFILE\n');
                        fprintf(fh,'%s\n',strcat('airfoils/',sect.airfoil));
                    end
                    fprintf(fh,'CLAF\n');
                    fprintf(fh,'%f\n',sect.CLAF);
                
                end
                
            end
            fclose(fh);
        end
        function []=avlplot(obj)
            obj.writeAVL('tempplane');
            fh=fopen('temprunfile.txt','w');
            fprintf(fh,'oper\n');
            fprintf(fh,'g\n');
            fprintf(fh,'tr\n');
            fprintf(fh,'k\n');
            fprintf(fh,'\n\n\nquit');
            fclose(fh);
            
            disp('Press spacebar to close AVL plot...');
            if ispc
                command=sprintf('avl tempplane.avl < temprunfile.txt > tempoutput.txt');
            elseif ismac
                command=sprintf('%s/avl3.35 tempplane.avl < temprunfile.txt > tempoutput.txt',pwd);
            end
            system(command);
        end
    end
end

