m_acc = 0.122;

mu_acier = 8000;
V_acier = 0.04^2*0.003;

A_colle = pi*0.02^2;
I_colle = pi*0.04^4/64;
h_colle = 0.002;
G_colle = 1e7; %http://www.tainstruments.com/pdf/literature/AAN001e_Hot_Melts.pdf

f1 = 1/(2*pi) * sqrt( G_colle*A_colle/h_colle / (m_acc + mu_acier*V_acier))

f2 = 1/(2*pi) * sqrt( G_colle*I_colle/h_colle / ((m_acc + mu_acier*V_acier) * 0.0466^2))