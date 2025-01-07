figure;
Rect = [0.19, 0.07, 0.775, 0.845];
% Rect = [0.1, 0.5, 0.7, 0.9];
AxisPos = moPlotPos(1, 4, Rect)
for i = 1:16
  axes('Position', AxisPos(i, :));
end