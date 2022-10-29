function out = detector(input)
symbols = ["00", "01", "10", "11"];
cons = [[sqrt(2)/2, sqrt(2)/2]
        [-sqrt(2)/2, sqrt(2)/2]
        [-sqrt(2)/2, -sqrt(2)/2]
        [sqrt(2)/2, -sqrt(2)/2]];
d = input - cons;
d = sqrt((d(:,1)).^2+(d(:,2)).^2);
[~, i] = min(d);
out = symbols(i);
end