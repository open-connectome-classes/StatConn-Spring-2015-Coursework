% akim1 150513

% grid
figure(1);
subplot(3,1,1);
imagesc(y_nonrand);
title('nonrand (spatial)');

subplot(3,1,2);
imagesc(y_rand);
title('rand (spatial)');

subplot(3,1,3);
imagesc(y_mix);
title('mix (spatial)');


% adjacency
figure(2);
subplot(3,2,1);
imagesc(mat_nonrand);
title('nonrand+orig (adjacency)');
subplot(3,2,2);
a = coarsen_conn(n, n, mat_nonrand);
imagesc(a);
title('nonrand+coarse (adjacency)');

subplot(3,2,3);
imagesc(mat_rand);
title('rand+orig (adjacency)');
subplot(3,2,4);
b = coarsen_conn(n, n, mat_rand);
imagesc(b);
title('rand+coarse (adjacency)');

subplot(3,2,5);
imagesc(mat_mix);
title('mix+orig (adjacency)');
subplot(3,2,6);
c = coarsen_conn(n, n, mat_mix);
imagesc(c);
title('mix+coarse (adjacency)');


% clustering
figure(3);
c_nonrand = reshape(kmeans(mat_nonrand, 4), n, n);
c_cnonrand = reshape(kmeans(a, 4), n/2, n/2);
c_rand = reshape(kmeans(mat_rand, 4), n, n);
c_crand = reshape(kmeans(b, 4), n/2, n/2);
c_mix = reshape(kmeans(mat_mix, 4), n, n);
c_cmix = reshape(kmeans(c, 4), n/2, n/2);

subplot(3,2,1);
imagesc(c_nonrand);
title('nonrand+orig (clustering)');
subplot(3,2,2);
imagesc(c_cnonrand);
title('nonrand+coarse (clustering)');

subplot(3,2,3);
imagesc(c_rand);
title('rand+orig (clustering)');
subplot(3,2,4);
imagesc(c_crand);
title('rand+coarse (clustering)');

subplot(3,2,5);
imagesc(c_mix);
title('mix+orig (clustering)');
subplot(3,2,6);
imagesc(c_cmix);
title('mix+coarse (clustering)');