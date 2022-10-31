% Total heat transfer rate of a rectangular enclosure
sigma = 5.6693e-8;
A = [3, 5, 3, 5]; epsilon = [0.7, 0.3, 0.85, 0.45];
T = [550, 700, 650, 600];
F = -[0, 0.3615, 0.277, 0.3615;...
0.2169, 0, 0.2169, 0.5662;...
0.277, 0.3615, 0, 0.3615;...
0.2169, 0.5662, 0.2169, 0];
Q = zeros(1, length(A));
c = zeros(1, length(A));
b = sigma*epsilon./(1-epsilon).*(1-c).*T.^4+c.*Q./A;
d = (1-c).*1./(1-epsilon)+c;
F = F+diag(d);
q0 = F\b';
Q = A.*epsilon./(1-epsilon).*(1-c).*(sigma*T.^4-q0')
q = Q./A