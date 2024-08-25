p_hole = 0.25;
x_hole = [p_hole]; %[mm]
y_hole = [p_hole]; %[mm]
r_hole = 20; %[mm]
w_hole = r_hole/10; %[mm] Beam waist of Gaussian
A = zeros(V,H);
for jhole = 1:length(x_hole)
    x_ap = x_hole(jhole);
    y_ap = y_hole(jhole);
    RHO = ((X-x_ap).^2 + (Y-y_ap).^2).^0.5;
    PHI = atan2((Y-y_ap),(X-x_ap));
    A(RHO <= r_hole) = exp(-(RHO(RHO <= r_hole)).^2/w_hole^2);
    
end
imagesc(A);
imshow(A);