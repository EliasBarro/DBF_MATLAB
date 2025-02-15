%Calculate T/W Needed for Cruise

function TW = cruise(WL, AR, CDmin, V, ALT)

q = dyn_P(ALT, V);

k = liftind_CD(AR);

TW = q*CDmin./WL + k./q*WL; %See Snorri 3-2

end