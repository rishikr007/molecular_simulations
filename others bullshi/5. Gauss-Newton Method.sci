clear;

x = 1:1:80;
height = 15;
centre = 50;
width = 10;

y = height * exp((-(x - centre).^2) / (2 * width^2));
y1 = y + 0.1 * (rand(1, 80) - 0.5) .* y; 

y1 = y1';

a = 10; b = 45; c = 7;

da = 0.1;  db = 0.1;  dc = 0.1; 

error_a = 1;  error_b = 1;  error_c = 1;

while (error_a > 0.0001 || error_b > 0.0001 || error_c > 0.0001)    
    Fx = a * exp((-(x - b).^2) / (2 * c^2));  
    fx = Fx';  

    Fxa = (a + da) * exp((-(x - b).^2) / (2 * c^2));
    Fxb = a * exp((-(x - (b + db)).^2) / (2 * c^2));
    Fxc = a * exp((-(x - b).^2) / (2 * (c + dc)^2));

    fxa = Fxa';  
    fxb = Fxb';
    fxc = Fxc';

    D = (y1 - fx);  

    z1 = (fxa - fx) / da;  
    z2 = (fxb - fx) / db;  
    z3 = (fxc - fx) / dc;  

    z = [z1 z2 z3];  

    lambda = 1e-6; 
    correction = (z' * z + lambda * eye(3)) \ (z' * D); 

    temp_a = a;
    temp_b = b;
    temp_c = c;

    a = a + correction(1);
    b = b + correction(2);
    c = c + correction(3);

    mprintf("Iter: a=%f, b=%f, c=%f\n", a, b, c);
    mprintf("Residual D:\n");
    disp(D);
    mprintf("Jacobian z:\n");
    disp(z);

    error_a = abs(a - temp_a) / abs(a);
    error_b = abs(b - temp_b) / abs(b);
    error_c = abs(c - temp_c) / abs(c);
end

mprintf("Final parameters: a =%f, b=%f, c=%f \n", a, b, c);

x_centered = x - mean(x);

plot(x_centered, y1, '*');    
yfit = a * exp((-(x - b).^2) / (2 * c^2));  
plot(x_centered, yfit, 'r');  
