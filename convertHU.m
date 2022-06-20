function RSP = convertHU(ct_orig,mvct,hu,y)
RSP=zeros(size(ct_orig));
for i=1:size(ct_orig,1)
    for j=1:size(ct_orig,2)
        for k=1:size(ct_orig,3)
            if mvct
                if ct_orig(i,j,k)<= -200
                    if ct_orig(i,j,k) > -1000
                        RSP(i,j,k) = (ct_orig(i,j,k) + 1000)*0.793/800;
                    else
                        RSP(i,j,k) = 0;
                    end
                elseif ct_orig(i,j,k)>140
                    RSP(i,j,k) = ct_orig(i,j,k)*2.281/4360 + 1.03;
                else
                    [~,ind] = find(hu==ct_orig(i,j,k));
                    RSP(i,j,k) = y(ind);
                end   
            else
                if ct_orig(i,j,k)<=1.1
                    RSP(i,j,k)=ct_orig(i,j,k)*0.0010457+1.0426;
                elseif 1.1<ct_orig(i,j,k)<=100
                    RSP(i,j,k)=ct_orig(i,j,k)*0.00037132+1.0893;
                elseif 100<ct_orig(i,j,k)
                    RSP(i,j,k)=ct_orig(i,j,k)*0.00055984+1.0302;
                end
            end
            
        end
    end
end
end