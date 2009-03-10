ii=1
while [  $ii -lt 121 ]; do
   source lancia.sh bs2phimumu_gen_${ii} bs2phimumu_reco_${ii}
   let ii=ii+1 
done


