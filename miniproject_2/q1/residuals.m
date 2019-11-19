function rls = residuals(x, t, y)
    rls = model(x,t)- y;
end