classdef avlsect
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        XYZle=[0 0 0];
        chord=1;
        ainc=0;
        Nspan=0;
        Sspace=0;
        airfoil='NACA0012';
        CLAF=1;
        
        controls=[];
        
    end
    
    methods
        function obj=avlsect(XYZle,chord,ainc,airfoil)
            obj.XYZle=XYZle;
            obj.chord=chord;
            obj.ainc=ainc;
            obj.airfoil=airfoil;
        end
    end
end

