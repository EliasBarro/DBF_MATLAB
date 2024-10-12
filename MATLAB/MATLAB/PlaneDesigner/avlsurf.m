classdef avlsurf
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name='tempsurface';
        Nchord=20;
        Cspace=1;
        Nspan=20;
        Sspace=1;
        Yduplicate=0; % what is the y-location of the axis of symmetry?
        angle=2; % incidence angle
        
        sections=[];
    end
    
    methods
        function obj=avlsurf(name,angle,sections)
            obj.name=name;
            obj.angle=angle;
            obj.sections=sections;
        end
    end
end

