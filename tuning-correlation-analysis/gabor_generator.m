%sizes = [7, 11, 14, 21];
%sizes = [7, 11, 14, 21, 35, 42, 49, 56, 60];

%orientations = [0 30 60 90 120 150];

%sizes=[20,30,40,50,60];
%orientations = 0:10:179;

sizes=[20,40,60]*2;
orientations = 0:12:180;

for size=sizes
    for o=orientations
        even = gabor(size, o, true);
        odd = gabor(size, o, false);
        
        %imwrite(even, sprintf('gabors/%d_%d_even.png', size, o));
        %imwrite(odd, sprintf('gabors/%d_%d_odd.png', size, o));

    end
end