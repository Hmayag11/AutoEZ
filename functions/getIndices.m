function ids0=getIndices(freq_mean,r)
id_delta=[];
id_theta=[];
id_alpha=[];
id_beta=[];
id_gamma=[];

for i=1:r
    if freq_mean(i)<4
        id_delta=[ id_delta i];
    elseif freq_mean(i)<8
        id_theta=[ id_theta i];
    elseif freq_mean(i)<12
        id_alpha=[id_alpha i];
    elseif freq_mean(i)<30
        id_beta=[ id_beta i];
    elseif freq_mean(i)<80
        id_gamma=[ id_gamma i];
    end
end

id_bb=1:r;
ids0{1}=id_delta;
ids0{2}=id_theta;
ids0{3}=id_alpha;
ids0{4}=id_beta;
ids0{5}=id_gamma;
ids0{6}=id_bb;
end
