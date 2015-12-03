
size(7cm,0);

pair[] n;
real h0 = 0;
real h1 = 0.5;
real h2 = 1;

// First tree
n[1] = (0,h0);
n[2] = (1,h0);
n[3] = (2,h0);
n[4] = (1.5,h1);
n[5] = (1,h2);

for (int i = 1; i <= 3; ++i) {
    dot(format("$%d$", i), n[i], S);
}
dot("$4$", n[4], NE);
dot("$5$", n[5], N);

draw(n[1]--n[5]);
draw(n[2]--n[4]);
draw(n[3]--n[5]);

// Second tree
n[1] = (3,h0);
n[2] = (4,h0);
n[3] = (5,h0);
n[4] = (4,h2);

for (int i = 1; i <= 3; ++i) {
    dot(format("$%d$", i), n[i], S);
}
dot("$4$", n[4], N);

draw(n[1]--n[4]);
draw(n[2]--n[4]);
draw(n[3]--n[4]);

// Third tree
n[1] = (6,h0);
n[2] = (7,h0);
n[3] = (8,h0);
n[4] = (6.5,h1);
n[5] = (7,h2);

for (int i = 1; i <= 3; ++i) {
    dot(format("$%d$", i), n[i], S);
}
dot("$4$", n[4], NW);
dot("$5$", n[5], N);

draw(n[1]--n[5]);
draw(n[2]--n[4]);
draw(n[3]--n[5]);
