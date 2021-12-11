Install supported module sources
==================================
Create folders for the modules you will use: `mkdir -p data/modules/`. If you plan to use several module sources, add subdirectories :
```bash
mkdir -p data/modules/BGSU
mkdir -p data/modules/RIN
mkdir -p data/modules/DESC
mkdir -p data/modules/JSON
mkdir -p data/modules/CSV
```

## CUSTOM JSON- OR CSV-FORMATTED MODULES
Just add you JSON-formatted modules to `data/modules/JSON/mydatabase.json`, according to the following format : 
```
{
    "1": {
        "sequence": "ACUAGCG&GGCUA&GU",
        "struct2d": "((((((.&.))))&))"
    },
    ...
}
```
You can use `'&'` to indicate sequence discontinuity, which leads to several components in the module.

You can also use CSV-formatted insertion sites (for example, obtained with Jar3d or BayesPairing) to `data/modules/CSV`, following one of these formats:

### The "BayesPairing" format:
Here k-loops can have any number of components k, you have to precise the start and end coordinates of each. The file should include the header.
```
Motif,Score,Start1,End1,Start2,End2...
motif1name,-19,29,38
motif2name,-28,71,80,90,96
...
```
Entries may not accumulate useless commas if they have a low number of components (don't `motif1name,-19,29,38,,`)

### The Jar3d format
Here the modules may only be 1-loops or 2-loops (HL or IL). There is a fixed number of columns per line, and undefined values are indicated with a dash `'-'`.
```
Motif,Rotation,Score,Start1,End1,Start2,End2
IL_43115.1,True,66,30,32,55,57
HL_35894.1,False,63,42,47,-,-
...
```

## CARNAVAL DATA (*Reinhartz et al, 2018*)

You first need to have the `unzip` command installed on your machine and the `networkx` package installed for Python 3. Then just run the script `Install_CaRNAval_RINs.py`.

If you have cloned the Git repository, just run :
```bash
cd scripts
python3 Install_CaRNAval_RINs.py
```
This will create files into `./data/modules/RIN/Subfiles`.

If not, or if you do not have the unzip command, download and extract manually the [CaRNAval dataset](http://carnaval.lri.fr/carnaval_dataset.zip) and place the files `RIN.py` and `CaRNAval_1_as_dictionnary.nxpickled` in the folder `data/modules/RIN/`, and run the python script.

*Note : CaRNAval is supposed to be a long-distance contact module dataset, not a SSE module dataset. It was supported for testing mostly, but you will not get the best performance from using it, it's not supposed to be loops.*

## THE RNA 3D MOTIF ATLAS DATA (*Petrov et al, 2013*, previously supported)
Source : see http://rna.bgsu.edu/rna3dhub/motifs/.

Get the latest version of the HL and IL module models from the [BGSU website](http://rna.bgsu.edu/data/jar3d/models/) and extract the Zip files. Put the HL and IL folders from inside the Zip files into `./data/modules/BGSU`. Note that only the latest Zip is required.

*Note : In Biorseo V1.0, you could use this modules directly because Biorseo was running Jar3d or BayesPairing for you. This is not the case anymore. You need to run these tools separately and get their results as a CSV file, see above how to format the CSV file.*

## RNA3DMOTIFS DATA (from the work of *Djelloul & Denise, 2008*, considered outdated)

If you use Rna3Dmotifs, you need to get RNA-MoIP's .DESC dataset: download it from [GitHub](https://github.com/McGill-CSB/RNAMoIP/blob/master/CATALOGUE.tgz). Put all the .desc from the `Non_Redundant_DESC` folder into `./data/modules/DESC`. Otherwise, you also can run Rna3Dmotifs' `catalog` program to get your own DESC modules collection from updated 3D data (download [Rna3Dmotifs](https://rna3dmotif.lri.fr/Rna3Dmotif.tgz)). You also need to move the final DESC files into `./data/modules/DESC`.