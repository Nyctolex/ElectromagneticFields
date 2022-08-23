%% Parameters
clear
clc
sections_num = 100;
txt = "";
id = [3 1 8 9 1 2 9 6 1];

[z1, z2, z3, z4, z5, z6, z7, z8, z9] = deal(id(1), id(2), id(3),  id(4), id(5), id(6), id(7), id(8), id(9));
V = 1;
d = 1+0.1*z3;
del = 0.01*(1+z1);
a = 5 + 0.5*z5;
e0 = 8.854e-12;
phi_total = [V*ones(sections_num,1); zeros(sections_num, 1)];

%% Question 2 section A1
%Plotting the plate capacitor via nomeric aproach.

x_len = 5*a;
m_point = 0;
x_start = m_point - x_len/2;
x_end = x_start + x_len;
x = linspace(x_start, x_end, sections_num);
y_surf = @(x) d*ones(1,length(x));
%A general function for the upper plate.
Z = calculate_Z_matrix(id, x_start, x_end, sections_num,V, y_surf);
%Getting the Z marix via a function.
heta_total = Z\phi_total;
%Solving the linear system.
heta_top = heta_total(1:end/2);
heta_buttom = heta_total(end/2+1:end);
heta_theo = e0*V/d;
%Theortial Heta of a plate capacitor.

figure;
hold on;
plot(x,heta_top);
plot(x,heta_buttom);
xlabel('x [m]');
ylabel('\eta [C/m^2]');
legend('Top Plate', 'Buttom Plate')

txt =print(txt, sprintf('Question 2 section A1: Theoretic Heta = e0*V/d =%s\n', double(heta_theo)));
image_name = sprintf('Q2 Section A1 - Surface charge density on plates (SN=%d)', sections_num);
suptitle(image_name);
hold off


%% Question 2 section A2
Q_total = trapz(x, heta_top);%sum(heta_top);%trapz(x, heta_top);
%Using the trapz function as a nomeric calculation of the integral.
capacity_nomeric = abs(Q_total)/V;
capacity_theo = e0*x_len/d;
%The theoretial (analitic) formula of the capacitance.

%printing all outputs:
txt =print(txt, sprintf('Question 2 section A2: Total Electric Charge (SN=%d): %s\n',sections_num, Q_total));
txt =print(txt, sprintf('Question 2 section A2: Capacity - Nomeric Calculation (SN=%d): %s\n', sections_num, capacity_nomeric));
txt =print(txt, sprintf('Question 2 section A2: Capacity - Analytic Calculation (SN=%d): %s\n', sections_num, capacity_theo));



%% Question 2 Section B1
%Calculating the electric dencity of the plates in nomeric aproach.

x_start = -a/4;
x_end = 5*a - a/4;
x = linspace(x_start, x_end, sections_num);
y_surf = @(x) d+del*cos(2*pi*x/a);

%A general function of the upper plate
Z = calculate_Z_matrix(id, x_start, x_end, sections_num,V, y_surf);
%Getting the Z marix via a function.
heta_total = Z\phi_total;
heta_top = heta_total(1:end/2);
heta_buttom = heta_total(end/2+1:end);
figure;
hold on;
plot(x,heta_top);
plot(x,heta_buttom);
xlabel('x [m]');
ylabel('\eta [C/m^2]');
image_name = sprintf('Q2 Section B1 - Surface charge density on plates (SN=%d)', sections_num);
suptitle(image_name);
hold off


%% Question 2 Section B2
%Verefing if the total electric charge of the plates is equal.
Q_top = trapz(x, heta_top);
Q_buttom = trapz(x, heta_buttom);
%Using the trapz function as a nomeric calculation of the integral.
txt = print(txt, sprintf('Question 2 section B2: Electric Charge of upper plane (SN=%d): %s\n', sections_num, Q_top));
txt = print(txt, sprintf('Question 2 section B2: Electric Charge of buttom plane (SN=%d): %s\n', sections_num, Q_buttom));
if (Q_top+Q_buttom) == 0
    txt = print(txt, sprintf('Question 2 section B2 conculution: The charge is equal and opposite\n'));
else
    txt = print(txt, sprintf('Question 2 section B2 conculution: The charge is *not* equal and opposite\n'));
    txt = print(txt, sprintf('Question 2 section B2 conculution: Diffrencr of %s\n',Q_top+Q_buttom));
    %Printing the differnce of charge between the plates.
end

%% Question 2 Section B3
%Calculating the capacity of the device via nomerical aproach.
Q_total = Q_top;
capacity_nomeric = abs(Q_total)/V;
txt = print(txt, sprintf('Question 2 section B3: Capacity - Nomeric Calculation (SN=%d): %s\n', sections_num, capacity_nomeric));


%% Question 2 Section C
resolution = 30;
%Number of b/d that we were asked to check.
syms A0 A1;
d_start = 2*del;
d_end = 10*del;
num_d = resolution;
d_scale = linspace(d_start, d_end, num_d);
capacity = zeros(length(d_scale), 1);
%A vector of the capacity via nomerical aprroach.
analytic_capacity = zeros(length(d_scale), 1);
%A vector of the analytic capacity as was calculated in question 1.
symbolic_capacity = Symbolic_Analytic_capcity(A0, A1,e0, a, V, x_start, x_end);
%A general formula for the analitical capacity.

for i = 1:length(d_scale)
    d = d_scale(i);
    y_surf = @(x) d+del*cos(2*pi*x/a);
    %General function of the upper plate.
    Z = calculate_Z_matrix(id, x_start, x_end, sections_num,V, y_surf);
    %Getting the Z marix via a function.
    heta_total = Z\phi_total;
    %Solving the linear system
    heta_top = heta_total(1:end/2);
    heta_buttom = heta_total(end/2+1:end);
    Q_top = trapz(x, heta_top);
    Q_total = Q_top;
    capacity_nomeric = abs(Q_total)/V;
    an = calculate_An_matrix(2, V, d, del, a);
    %getting the Augmentes of the analytic series.
    capacity(i) = capacity_nomeric ;
    analytic_capacity (i) = subs(symbolic_capacity, {A0 A1}, {an(1), an(2)});
end
section_c_err = abs(capacity - analytic_capacity)./capacity;
%calculating the relative error between the nomeric and analytic capacity
%as a function of d
txt = print(txt, sprintf('Question 2 section E: Error between analytic and nomerical calculation of section C (SN=%d):\n', sections_num));
txt = print_vector(txt, section_c_err);
txt = print(txt, sprintf('Question 2 section E: Avrage error between analytic and nomerical calculation of section C (SN=%d)%s\n',sections_num, mean(section_c_err)));
figure;
hold on;
plot(d_scale(2:end), capacity(2:end));
plot(d_scale(2:end), analytic_capacity(2:end));
legend('nomerical', 'analytic');
xlabel('d [m]');
ylabel('capacity[F]');
image_name = sprintf('Q2 Section C - Capacity as function of d (SN=%d)', sections_num);
suptitle(image_name);
hold off


%% Question 2 Section D
b_num = resolution;
b_scale = linspace(0, 1, b_num);
d = 1+0.1*z3;
y_surf = @(x) d+del*cos(2*pi*x/a);

an = calculate_An_matrix(2, V, d, del, a);
%getting the Augmentes of the analytic series.
capacity = zeros(length(b_scale), 1);
%A vector of the capacity via nomerical aprroach.
analytic_capacity = zeros(length(b_scale), 1);
for i = 1:length(b_scale)
    b = b_scale(i);
    x_start =  b*a/2;
    x_end = 5*a+b*a/2;
    Z = calculate_Z_matrix(id, x_start, x_end, sections_num,V, y_surf);
    %Getting the Z marix via a function.
    heta_total = Z\phi_total;
    %Solving the linear system
    heta_top = heta_total(1:end/2);
    heta_buttom = heta_total(end/2+1:end);
    Q_top = trapz(x, heta_top);
    Q_buttom = trapz(x, heta_buttom);
    %Using the trapz function as a nomeric calculation of the integral.
    Q_total = Q_top;
    capacity_nomeric = abs(Q_total)/V;
    capacity(i) = capacity_nomeric ;
    symbolic_capacity = Symbolic_Analytic_capcity(A0, A1,e0, a, V, x_start, x_end);
    analytic_capacity (i) =subs(symbolic_capacity, {A0 A1}, {an(1), an(2)});
end
section_d_err = abs(capacity - analytic_capacity)./capacity;

txt = print(txt, sprintf('Question 2 section E: Error between analytic and nomerical calculation of section D (SN=%d):\n', sections_num));
txt = print_vector(txt, section_d_err);
%calculating the relative error between the nomeric and analytic capacity
%as a function of b
txt = print(txt, sprintf('Question 2 section E: Avrage error between analytic and nomerical calculation of section D (SN=%d): %s\n',sections_num, mean(section_d_err)));

%plotting
figure;
hold on;
plot(b_scale, capacity);
xlabel('b');
ylabel('capacity[F]');
image_name = sprintf('Q2 Section D - Nomeric Capacity as function of b (SN=%d)', sections_num);
suptitle(image_name);
hold off



figure;
hold on;
plot(b_scale, capacity);
plot(b_scale, analytic_capacity);
legend('nomerical', 'analytic');
xlabel('b');
ylabel('capacity[F]');
image_name = sprintf('Q2 Section D - Capacity as function of b (SN=%d)', sections_num);
suptitle(image_name);
hold off



%% Question 2 Section E
%plot section C error
figure;
hold on;
plot(d_scale, section_c_err);
xlabel('d');
ylabel('Relative error');
image_name = sprintf('Q2 Section E - relative error of capacity as function of d  (SN=%d)', sections_num);
suptitle(image_name);
hold off

%plot section D error
figure;
hold on;
plot(b_scale, section_d_err);
xlabel('b');
ylabel('Relative error');
image_name = sprintf('Q2 Section E - relative error of capacity as function of b (SN=%d)', sections_num);
suptitle(image_name);
hold off


%% functions
function [dm] = calculate_dm(y_surf, x_start, x_end, sections_num)
%Calculating delta_m as a function of the length of the plate (y_suf)
%divided by the number of section number we have chosen.
    syms t;
    dl = sqrt(diff(t, t)^2 + diff(y_surf(t), t)^2);
    L = int(dl, t, x_start, x_end);
    dm = double(L / sections_num) ;
end

function [Z] = calculate_Z_matrix(id, x_start, x_end, sections_num,V, y_surf)
%The function is buildiung the Z matrix and returns it.
    [z1, z2, z3, z4, z5, z6, z7, z8, z9] = deal(id(1), id(2), id(3),  id(4), id(5), id(6), id(7), id(8), id(9));
    am = 1e6;
    e0 = 8.854e-12;
    dm = calculate_dm(y_surf, x_start, x_end, sections_num);
    x = linspace(x_start, x_end, sections_num);
    y_total = [y_surf(x), zeros(1,sections_num)];
    circular_indexing = @(index, len) (1 + mod(index-1, len));
    %since we are using one matrix for the linear system but two plates,
    %circular_indexing return the matching index when the x index is out
    %of range.
    r = @(n) [x(circular_indexing(n, sections_num)), y_total(n)];
    zmn = @(m, n) (-(dm/(2*pi*e0))*log((rssq(r(n)-r(m))+double(m==n))/am));
    zmm = -(dm/(2*pi*e0))*(log(dm/(2*am)) -1);
    z_func = @(m,n) zmn(m,n)*(m~=n)+zmm*(m==n);
    %chooses between the diffrent formulas of a Z cell (if m is equal or not to n).
    Z = zeros(sections_num*2 , sections_num*2);
    for m=1:2*sections_num
        for n = 1:2*sections_num
            Z(m, n) = z_func(m,n);
        end
    end
end

function [capacity_theo] = Symbolic_Analytic_capcity(A0, A1,e0, a, V, x_start, x_end)
%Calculating the theoritical capacity according the Question section F in
%symbolic form and according to Boundary value problem
    syms x y;
    phi = A0*y+A1*cos(2*pi*x/a)*sinh(2*pi*y/a);
    E1 = [-diff(phi, x), -diff(phi, y)];
    E1 = subs(E1, y, 0);
    E2 = [0, 0];
    heta = e0*(E1 - E2);
    heta = heta(2);
    Q_total = int(heta, x, x_start, x_end);
    capacity_theo = abs(Q_total)/V;
end


function [an] = calculate_An_matrix(N, V, d, del, a)
%A function which calculate the An series from index 1 to N.
%Notice that this is the same function as in Question 1
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

function [txt] = print_vector(txt, object)
%Displating the vector and adding its string into txt.
    temp_txt = sprintf("\n%s\n", object);
    txt = txt+temp_txt;
    disp(object);
end
