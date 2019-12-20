# mimt
MIami Multi Track tools

save_geodmod.py generates data that used by Geodmod software. "save_geodmod.py --help" can give the detail message about input parameters you have to provide. 

save_gbis_mimt.py generates data that used by GBIS software. "save_gbis_mimt.py --help" can give the detail message about input parameters you have to provide.

generate_horzvert.py generates horizontal and vertical data file. "generate_horzvert.py --help" can give the detail message about input parameters you have to provide.

save_for_modelling.py uses template file to generate data used by Geodmod/GBIS or generate horizontal/vertical files. The example template file is Darbandikan.txt. Help massage give the detail about the usage.

# Access to mintpy data products
Browse through all directories using: https://129.114.104.223/data/HDF5EOS 

Get them via wget using:
```
wget -r -nH --cut-dirs=2 --no-parent --reject="index.html*" --no-check-certificate https://129.114.104.223/data/HDF5EOS
wget -r -nH --cut-dirs=2 --no-parent --reject="index.html*" --no-check-certificate https://129.114.104.223/data/HDF5EOS/KashgarSenAT129

```
