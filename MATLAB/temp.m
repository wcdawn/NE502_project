z = linspace(0,H);
qpp = zeros(length(z),1);
for i = 1:length(z)
    qpp(i) = q0pp * ((pi * (H + lambda - z(i))) / (He)) * sin((pi * (H + lambda - z(i))) / (He));
end