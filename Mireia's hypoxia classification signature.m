% Calculate TBRmax
structNumV = [1];
badSlicesStr = "None";
[planC,hiSuvCtStrV,hiSuvPetStrV] = ...
   preprocessMAASTROscans(structNumV,badSlicesStr,planC);
featureS = getHypoxiaCrispinOrtuzarFeatures(hiSuvCtStrV,hiSuvPetStrV,planC);
TBRmax = NaN*ones(1,length(hiSuvCtStrV));
for i=1:length(hiSuvCtStrV)
  if getStructureVol(hiSuvCtStrV(i),planC) > 10
    P90_norm = (featureS(i).P90-11.5446)/5.1479;
    LRHGLE_norm = (featureS(i).lrhgle - 7254.5)/1122.1;
    TBRmax(i) = 1.9061 + 0.32381*P90_norm+0.13032*LRHGLE_norm;
  end
end
