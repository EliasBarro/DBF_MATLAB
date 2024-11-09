%Calculate Lift Induced Drag Constant

function k = liftind_CD(AR)

e = spn_Eff(AR); %Calculate Oswald's Span Efficiency

k = 1/(3.14*AR*e); %See Snorri Eqn. 15-7

end