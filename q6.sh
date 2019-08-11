. ~/.bashrc
unset SGE_ROOT

for d in */;
do
# cp curie.f "$d";
# cp erg.f "$d";
# cp temperature.f "$d";
cp mani6.py "$d";
cp q6.sh "$d";
cp rr6.sh "$d";
done


 

for subdir in */;
do 
cd "$subdir"
bsub < rr6.sh
cd ..
done


