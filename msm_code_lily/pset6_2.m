% Problem 2

%% 2a)
load clust_data.dat

[IDX, center] = kmeans(clust_data,2);
center

%% 2b)
clust1 = zeros(10,2);
clust2 = zeros(10,2);
n = 1;
m = 1;
hold on
for i=1:200
    if IDX(i)==1
        clust1(n,:) = clust_data(i,:);
        n = n+1;
    else
        clust2(m,:) = clust_data(i,:);
        m = m+1;
    end
end

%% 2c)
plot(clust1(:,1),clust1(:,2), "b.")
plot(clust2(:,1),clust2(:,2), "g.")

hold on
plot(center(1,1),center(1,2),"ko")
plot(center(2,1), center(2,2), "ko")

legend("cluster 1, cluster 2, center 1, center 2")
xlabel("X")
ylabel("Y")