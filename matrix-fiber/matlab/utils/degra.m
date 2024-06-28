syms d p a1 a2

wd = (1 - d)^p / ((1 - d)^p + a1 * d + a1 * a2 * d^2);

wd_prime = diff(wd, d);

% 显示结果
% pretty(wd);
% pretty(wd_prime);

% 数值替换
p_value = 2;
a1_value = 80;
a2_value = -0.5;

% 定义 wd 的匿名函数
wd_fun = matlabFunction(subs(wd, [p, a1, a2], [p_value, a1_value, a2_value]));
wd_d_fun = matlabFunction(subs(wd_prime, [p, a1, a2], [p_value, a1_value, a2_value]));
% 绘制 wd 的图像
figure;
fplot(wd_fun, [0, 0.5]);
figure();
fplot(wd_d_fun, [0, 0.5]);
xlabel('d');
ylabel('wd');
title('Plot of wd');
grid on;


