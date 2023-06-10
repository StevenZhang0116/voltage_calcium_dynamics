r1 = load('meanfield_result_2002.mat').result;
r2 = load('meanfield_result_2011_0.mat').result;

diff1 = max(r1(:,3)) - min(r1(:,3));
diff2 = max(r2(:,3)) - min(r2(:,3));

s1 = r1(:,8);
s2 = r2(:,8);

area1 = 0; area2 = 0;
for i = 1:length(s1)
    area1 = area1+s1(i);
end
area1 = area1 + s1(1) + s1(end);
for i = 1:length(s2)
    area2 = area2+s2(i);
end
area2 = area2 + s2(1) + s2(end);

figure()
hold on
plot(r1(:,1),r1(:,8))
plot(r2(:,1),r2(:,8))
title("SR Ca^{2+} Release Flux")
xlabel("Time(ms)")
ylabel("SR Ca^{2+} Release Flux")
lgd = legend('Adjusted 2002','2011');
lgd.FontSize = 14;
hold off