size(11cm,0);
import geometry;

pair[] n;
real h0 = 0;
real h1 = 0.5;
real h2 = 1;
real h3 = 2.5;
real hl = -1;

triangle t = rotate(225) * triangle(b=0.33,alpha=90,c=0.33);

// Nodes are numbered left to right, root downwards
n[1] = (2.5, h3);
n[2] = (1, h2);
n[3] = (4, h2);
for (int i = 0; i < 5; ++i) {
    if (i < 4) {
        n[i + 4] = (i, h0);
    } else {
        n[i + 4] = (i + 1, h0);
    }
    draw(shift(n[i + 4]) * t);
}
n[9] = (0.5, h1);
draw(n[2]--n[1]--n[3], dashed);
draw(n[4]--n[2]--n[6]);
draw(n[7]--n[3]--n[8]);
draw(n[5]--n[9]);

dot("$z$", n[1], N);
dot("$x$", n[2], NW);
dot("$y$", n[3], NE);
dot("$a$", n[4], NW);
dot("$b$", n[5], N);
dot("$c$", n[6], N);
dot("$d$", n[7], N);
dot("$e$", n[8], N);
dot("$\alpha$", n[9], NW);

label("(i)", (2.5, hl));

// The intermediate state.

int base = 7;
real off = 0.75;
draw((base - off, 0)--(base - off, h3), gray);

n[1] = (base + 2.5, h3);
n[2] = (base + 1, h2);
n[3] = (base + 4, h2);
for (int i = 0; i < 5; ++i) {
    n[i + 4] = (base + i + 0.5, h0);
    draw(shift(n[i + 4]) * t);
}
draw(n[2]--n[1]--n[3], dashed);

dot("$z$", n[1], N);
dot("$x$", n[2], NW);
dot("$y$", n[3], NE);
dot("$a$", n[4], N);
dot("$b$", n[5], N);
dot("$c$", n[6], N);
dot("$d$", n[7], N);
dot("$e$", n[8], N);

label("(ii)", (base + 2.5, hl));

// The final state.
base = 14;
real off = 1.25;
draw((base - off, 0)--(base - off, h3), gray);

n[1] = (base + 2.5, h3);
n[2] = (base + 1, h2);
n[3] = (base + 4, h2);
for (int i = 0; i < 5; ++i) {
    if (i == 1) {
        n[i + 4] = (base + 4, h0);
    } else if (i < 4) {
        n[i + 4] = (base + i, h0);
    } else if (i == 4) {
        n[i + 4] = (base + i + 1, h0);
    }
    draw(shift(n[i + 4]) * t);
}
n[9] = (base + 4.5, h1);
draw(n[2]--n[1]--n[3], dashed);
draw(n[4]--n[2]--n[6]);
draw(n[7]--n[3]--n[8]);
draw(n[5]--n[9]);

dot("$z$", n[1], N);
dot("$x$", n[2], NW);
dot("$y$", n[3], NE);
dot("$a$", n[4], N);
dot("$b$", n[5], N);
dot("$c$", n[6], N);
dot("$d$", n[7], N);
dot("$e$", n[8], NE);
dot("$\beta$", n[9], NE);

label("(iii)", (base + 2.5, hl));

