function eps=die(modee,freq)
e=1.60217662*10^-19;
eps0=8.85*10^-12;
b=modee*e;
vol=(8.386494280/2*10^-10)^3;
f=2.9980*10^8*freq*100;
eps=b'*b/4/pi^2/eps0/f^2/vol*6.02*10^(23)/(10^-3);
end
