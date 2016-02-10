function  drawResult(input_name)

control_name = strcat(input_name, '_controls.txt');
spline_name = strcat(input_name, '_spline.txt');

P = load(input_name);
figure;
plot(P(:,1), P(:,2),'b.','MarkerSize',15);
%plot(P(:,1), P(:,2),'b');

hold on;
C = load(control_name);
C2 = [C; C(1,:)];
plot(C2(:,1),C2(:,2),'r','LineWidth', 1);
plot(C2(:,1),C2(:,2), 'rx', 'MarkerSize',10);


S = load(spline_name);
S2 = [S; S(1,:)];
plot(S2(:,1),S2(:,2),'g','LineWidth', 2);
end








