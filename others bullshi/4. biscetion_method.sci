function y = f(x)
    y = x^2 - 4; 
endfunction

function root = bisection(f, a, b, tol, max_iter)
    if f(a) * f(b) >= 0 then
        error("The function must have different signs at the endpoints a and b.");
        return;
    end

    for iter = 1:max_iter
        c = (a + b) / 2;
        
        if abs(f(c)) < tol then
            break;
        elseif f(a) * f(c) < 0 then
            b = c; 
        else
            a = c;
        end
    end
    root = c;
    disp("Iterations: " + string(iter));
endfunction

a = input("Lower limit: ");
b = input("Upper limit: ");
tol = input("Tolerance: ");
max_iter = 100;

root = bisection(f, a, b, tol, max_iter);

disp("The root is: " + string(root));
