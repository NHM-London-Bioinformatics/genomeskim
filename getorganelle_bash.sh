

# Activate environment
conda activate getorganelle

# Set the database required
Fdb="animal_mt"

# Set the default databases location (conda only)
dbpath=$(conda info --envs | grep "\*" | sed -e "s/^.* //")
dbpath=$(echo "$dbpath/getorganelledb/")
mkdir -p "$dbpath"
export GETORG_PATH="$dbpath"

# getorganelle command
cmd="get_organelle_from_reads.py -1 ../../testdata/minimal/ptilo.R1.fq.gz -2 ../../testdata/minimal/ptilo.R2.fq.gz -F animal_mt -o getorganelle_ptilo_test/"


# TODO: run getorganelle - if fails, capture output and parse for missing databases

# Install missing databases
if [[ $(get_organelle_config.py --list | wc -l) != 14  ]]
then
   #get_organelle_config.py -a all
   # This currently fails because the sha256sums are wrong. Use this 
   git clone https://github.com/Kinggerm/GetOrganelleDB
   get_organelle_config.py -a all --use-local ./GetOrganelleDB/0.0.1/
   # TODO write a more sophisticated error catch for these two options
   
   # TODO then re-run getorganelle
fi




