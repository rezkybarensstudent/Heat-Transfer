% Transient Heat Conduction Time and Temperature Distributions in a Semi-Infinite Solid

% Range of the temperature in the semi-infinite solid (1)
tau = linspace(0.01, 3, 30);
% Range of the temperature in the semi-infinite solid (2)
eta = linspace(0, 5, 20);
% Generate x and t matrices for three-dimensional plots
[x, t] = meshgrid(eta, tau);
%{
The inline command lets you create a function of any number of variables
by giving a string containing the function followed by a series of strings
denoting the order of the input variables.
%}
theta = inline('erfc(0.5*x./t)-exp(x+t.^2).*erfc(0.5*x./t+t) ', 'x', 't');
% Create a figure graphics object
figure(1)
% Mesh surface plot
mesh(x, t, theta(x, t))
% eta Label
xlabel('\eta')
% tau Label
ylabel('\tau')
% theta Label
zlabel('\theta')
% Create a figure graphics object (2)
figure(2)
% eta Number
eta = 0:5;
% Range tau Number
tau = linspace(0.01, 4, 40);
for k = 1:length(eta)
thet = theta(eta(k), tau);
plot(tau, thet, 'k-')
text(.92*4,1.02*thet(end), ['\eta = ' num2str(eta(k))])
hold on
end
% tau Label (2)
xlabel('\tau')
% theta Label (2)
ylabel('\theta')
