function y = energy_to_peakPos(energy,data,str)

data_energy = zeros(1, length(data));
for i = 1 : length(data)
    data_energy(i) = data(i).energy;
end

[~,ind] = min(abs(data_energy-energy));

if strcmp(str,'range')
    y = data(ind).range;
else
    y = data(ind).peakPos;
end

end