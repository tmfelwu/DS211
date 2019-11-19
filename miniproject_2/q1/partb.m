load('data.mat')
tic()
sol = gn(@model, @jacobian_f, @residuals, @jacobian, t, y, 1, 200);
toc()
sol