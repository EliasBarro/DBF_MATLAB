% DO NOT use the integer values these enumerations are assigned to, that
% % defeats the purpose and makes things less clear. The integers are just
% there so that you can plot mission phase as a timehistory
classdef mp < uint32
    enumeration
        terminate   (42)   % - 
        goto        ( 0)   % phase to go to
        config      (47)   % statements to run
        
        % ground modes
        ground_thr  ( 1)   % thr setting (0-1)
        ground_rot  ( 2)   % rotate_time
        
        % flight modes
        to_trans    ( 3)   % time to hold (sec), hold last thr
        cruise_thr  ( 4)   % thr setting (0-1)
		hold_spd    ( 5)   % time to hold (sec)
        climb_max   ( 6)   % gamma (rad)             - max thr
        glide       ( 7)   % gamma (rad)             -  no thr
        opt_climb   ( 8)   % CL                      - max thr
        opt_glide   ( 9)   % CL                      -  no thr
        turn_gload  (10)   % g-load of turn, hold last thr
        turn_radi   (11)   % radius of turn (m), hold last thr

        cruise_thr_down(12)
        
    end
    methods
        function mp_str = get_str(obj)
			if     obj == mp.terminate  , mp_str = 'terminate' ; %#ok<ALIGN>
			elseif obj == mp.goto       , mp_str = 'goto'      ;
            elseif obj == mp.config     , mp_str = 'config'    ;
			elseif obj == mp.cruise_thr , mp_str = 'cruise_thr';
			elseif obj == mp.hold_spd   , mp_str = 'hold_spd'  ;
			elseif obj == mp.turn_gload , mp_str = 'turn_gload';
			elseif obj == mp.turn_radi  , mp_str = 'turn_radi' ;
			elseif obj == mp.climb_max  , mp_str = 'climb_max' ;
			elseif obj == mp.glide      , mp_str = 'glide'     ;
			elseif obj == mp.opt_climb  , mp_str = 'opt_climb' ; % not valid yet!
			elseif obj == mp.opt_glide  , mp_str = 'opt_glide' ; % not valid yet!
			elseif obj == mp.ground_thr , mp_str = 'ground_thr';
            elseif obj == mp.ground_rot , mp_str = 'ground_rot';
            elseif obj == mp.to_trans   , mp_str = 'to_trans'  ;
            elseif obj == mp.cruise_thr_down, mp_str = 'cruise_thr_down';
            else
                disp('wut, probs were had');
            end     
        end
    end
end

% one day should condense cruisethr/climb/glide combos into fly gamma and fly CL 
% where thr and gamma/CL are TWO inputs in a cell array

