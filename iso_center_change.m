function center = iso_center_change(iso_center,orig_nii)

center = zeros(size(iso_center));
for i = 1 : size(iso_center,1)
    center(i,3) = int16((iso_center(i,3) - orig_nii.raw.srow_z(4)) / orig_nii.raw.srow_z(3));
    center(i,2) = int16((iso_center(i,2) - orig_nii.raw.srow_y(4)) / orig_nii.raw.srow_y(2));
    center(i,1) = int16((iso_center(i,1) - orig_nii.raw.srow_x(4)) / orig_nii.raw.srow_x(1));
end

end