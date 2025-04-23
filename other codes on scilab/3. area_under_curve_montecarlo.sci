a = input("Enter lower limit: ");
b = input("Enter upper limit: ");
n = input("Enter number of intervals: ");

deff("y = f(x)", "y = sin(x)");

n_hits = 0;
rect_area = (b - a) * 1; 

for i = 1:n
    x = a + rand() * (b - a);   
    y = rand();            
    if y < f(x) then
        n_hits = n_hits + 1;
    end
end

area = rect_area * n_hits / n;
disp(area);
