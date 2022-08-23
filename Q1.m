%% Parameters
clc
clear
close all

max_N_assumption = 4;
%The maximum N that for N+1 we would assume a_(n+1)=0
txt_path = sprintf("/output.txt");
%path for the txt output file
txt = "";
id = [3 1 8 9 1 2 9 6 1];
[z1, z2, z3, z4, z5, z6, z7, z8, z9] = deal(id(1), id(2), id(3),  id(4), id(5), id(6), id(7), id(8), id(9));

V = 1;
d = 1+0.1*z3;
del = 0.01*(1+z1);
a = 5 + 0.5*z5;
e0 = 8.854e-12;
%% Question 1 Section D1
%Draw potrntial on upper surf
N=max_N_assumption;
an = calculate_An_matrix(N, V, d, del, a);
figure;
x = linspace(0, 5*a, 300);
y = d+del*cos(2*pi*x/a);
sum_series = 0;
%the sum of the protential formula except for An[0]*y
for j = 1:N-1
    sum_series =sum_series+ an(j+1)*cos(2*pi*j.*x/a).*sinh(2*pi*j.*y/a);
end
potential = an(1).*y+sum_series;

hold on;
plot(x, potential);
ylabel("Potential[V]");
xlabel("x[m]");
potential_error = max(potential) - 1;
%displaying the nomerical answers ans saving it to the output.txt file:
txt = print(txt, sprintf('Question 1 section D1 An Vector:\n'));
txt = print_vector(txt, an);
txt = print(txt, sprintf('Question 1 section D1 conculution: Potential Error: %d\n',potential_error));
title = sprintf("Q1 section D1 - Potential on upper board (N=%d)", N);

suptitle(title);


%% Question 1 Section D2
N=11;
an = calculate_An_matrix(N, V, d, del, a);
div_an = @(n)an(n+1)/an(n);
%general formula of divination
quotient = [];
for i = 1:N-1
    quotient(i) = div_an(i);
end
n_scale = 0:1:N-2;
figure;
hold on;
plot(n_scale, abs(quotient));
plot(n_scale, ones(1, length(quotient))*exp(-2*pi*d/a));
xlabel("n");
ylabel("a_n_+_1/a_n")
image_name = sprintf('Question1 section D1 - Decay of coefficient as a function of n (N=%d)', N);
suptitle(image_name);
txt = print(txt, sprintf('Question 1 section D2: quotient vector of an+1/an:\n'));
txt = print_vector(txt, quotient);

%% Question 1 Section E
%calculating the potential in space.
%Initiation the parameter for this section
N= max_N_assumption;
an = calculate_An_matrix(N, V, d, del, a);

x = linspace(0, 5*a);
y = transpose(linspace(0, d+del));
sum_series = 0;
%the sum of the protential formula except for An[0]*y
for j = 1:N-1
    sum_series =sum_series+ an(j+1)*cos(2*pi*j.*x/a).*sinh(2*pi*j'*y/a);
end
potential = an(1).*y+sum_series;
figure;
contourf(x, y, potential, 30);
hold on;
plot(x, d+del*cos(2*pi*x/a),'b--', 'linewidth',3)
colormap(jet(100));
cb = colorbar;

set(gca,'YDir','normal');
ylabel("y[m]");
xlabel("x[m]");
image_name = sprintf("Question1 section E - Potential in space for (N=%d)", N);
suptitle(image_name);

%% Question 1 - section F
%Calculating a formuala for the capacity with symbols.

%Initiation the values to substatue.
subs_a = 5+z5/2;
subs_e0 = 8.854e-12;
subs_V = 1;
subs_d = 1+0.1*z3;
%Defning the symbols
syms x y A0 A1 e0 V a d A;
phi = A0*y+A1*cos(2*pi*x/a)*sinh(2*pi*y/a);
E1 = [-diff(phi, x), -diff(phi, y)];
E1 = subs(E1, y, 0);
E2 = [0, 0];
heta = e0*(E1 - E2);
%Boundary value problem
heta = heta(2);
%The orthogonal vector is y.


txt = print(txt, sprintf('Question 1 section F: Heta Formula:\n'));
txt = print_object(txt, heta);

%substituting paramaters.
x_start = 0;
x_end = 5*a;
Q_total = int(heta, x, x_start, x_end);
capacity_theo = Q_total/V;
capacity_subs = subs(capacity_theo, {A0 e0 V a}, {an(1) subs_e0 subs_V subs_a});
plate_capacitor = e0*A/d;

%Dispalying the calculation output and saving it to output.txt file
txt = print(txt, sprintf('Question 1 section F: Capacity Formula:%s \n', capacity_theo));
txt = print(txt, sprintf('Question 1 section F: Capacity=%s \n', double(capacity_subs)));
plate_capacitor_subs = subs(plate_capacitor, {A e0 d}, {5*subs_a subs_e0 subs_d});
txt = print(txt, sprintf('Question 1 section F: Plate Capacity Formula:%s \n', plate_capacitor));
txt = print(txt, sprintf('Question 1 section F: Plate Capacity=%s \n', double(plate_capacitor_subs)));


%% Question 1 - section G
d = 1+0.1*z3;
del = d/10;
V = 1;
a = 5 + 0.5*z5;
N=100;
an = calculate_An_matrix(N, V, d, del, a);
x = linspace(0, 5*a);
y = transpose(linspace(0, d+del));
sum_series = 0;
%the sum of the protential formula except for An[0]*y
for j = 1:N-1
    sum_series =sum_series+ an(j+1)*cos(2*pi*j.*x/a).*sinh(2*pi*j'*y/a);
end
potential = an(1).*y+sum_series;
figure;
contourf(x, y, potential, 30);
hold on;
plot(x, d+del*cos(2*pi*x/a),'b--', 'linewidth',3)
colormap(jet(100));
cb = colorbar;

set(gca,'YDir','normal');
ylabel("y[m]");
xlabel("x[m]");
image_name = sprintf("Question1 section G - Potential in space for (N=%d)", N);
suptitle(image_name);

%% functions
function [an] = calculate_An_matrix(N, V, d, del, a)
%A function which calculate the An series from index 1 to N.
    AN_matrix = sym(zeros(N,N));
    %The Augmented matrix of the linear system equation.
    an_minus = @(n)(pi*del/a)*(n-1)*cosh((2*pi*d/a)*(n-1));
    %a_(n-1) formula
    an = @(n)sinh(2*pi*n*d/a);
    %an formual
    an_plus = @(n)(pi*del/a)*(n+1)*cosh(2*pi*d*(n+1)/a);
    %a_(n+1) formula
    
    %first row
    AN_matrix(1, 1) = d;
    AN_matrix(1, 2) = (del*pi)/a*cosh(2*pi*d/a);
    %second row
    AN_matrix(2, 1) = del;
    AN_matrix(2, 2) = sinh(2*pi*d/a);
    AN_matrix(2, 3) = (2*pi*del/a)*cosh(4*pi*d/a);
    for j = 3:N
        AN_matrix(j, j-1) = an_minus(j-1);
        AN_matrix(j, j) = an(j-1);
        AN_matrix(j, j+1) = an_plus(j-1);
    end
    AN_matrix(N+1, N+1) = 1;
    AN_matrix(N, N+1) = 0;
    b = [[V]; zeros(N, 1)]; %answer vector of the linear system A*x = b
    an = AN_matrix\b;
    an = an(1:end-1);
    an = double(an);
end


function [txt] = print(txt, temp_txt)
%displaying temp_txt and adding it to txt
    txt = txt+temp_txt;
    disp(temp_txt);
end

function [txt] = print_object(txt, object)
%displaying object in 'pretty' form and adding its string value to txt
    temp_txt = sprintf("\n%s\n", object);
    txt = txt+temp_txt;
    pretty(object);
end

function [txt] = print_vector(txt, object)
%Displating the vector and adding its string into txt.
    temp_txt = sprintf("\n%s\n", object);
    txt = txt+temp_txt;
    disp(object);
end

