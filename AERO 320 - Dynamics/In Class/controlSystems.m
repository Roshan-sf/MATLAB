%gp = tf([1.515,0.1774],[1,0.739,0.921,0]);
gp = tf(1,[10,0,0]);

rlocus(gp)
%controlSystemDesigner(gp)