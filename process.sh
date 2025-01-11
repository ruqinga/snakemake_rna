
geofetch -i GSE162584 -n metadata --just-metadata
nohup geofetch -i GSE162584 > output.log &
mkdir sra
mv SRR* sra/

