deff('y=f(x)', 'y=sin(x)');

function integ = Trape(a, b, n, f)
    h = (b - a) / n;
    integ = (f(a) + f(b)) / 2;
    for i = 1:n - 1
        x = a + (i * h);
        integ = integ + f(x);
    end
    integ = integ * h;
endfunction

a = input("Enter lower limit: ");
b = input("Enter upper limit: ");
n = input("Enter number of intervals: ");

//disp(Trape(a, b, n, f));

if n <= 0 then
    disp("Error: The number of intervals must be a positive integer.");
else
    result = Trape(a, b, n, f);
    printf("The approximate value of the integral from %f to %f is: %f\n", a, b, result);
end
