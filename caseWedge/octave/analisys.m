
x = csvread("./x.csv");
y = csvread("./y.csv");

r = csvread("./r.csv");
ru = csvread("./ru.csv");
rv = csvread("./rv.csv");
rE = csvread("./rE.csv");

%colormap ("default");


contourf(rv);
axis("equal")


