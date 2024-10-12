%Calculate Oswald's Span Efficiency

function e = spn_Eff(AR)

e = 1.78*(1-0.045*AR^0.68)-0.64; %Raymer's estimation for straight wings; see Snorri 9.5.14

end