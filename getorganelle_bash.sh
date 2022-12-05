

# Activate environment
# conda create --name getorganelle -c bioconda getorganelle
conda activate getorganelle

# Set the default databases location (conda only)
dbpath=$(conda info --envs | grep "\*" | sed -e "s/^.* //")
dbpath=$(echo "$dbpath/getorganelledb/")
mkdir -p "$dbpath"
export GETORG_PATH="$dbpath"

# getorganelle command
cmd="get_organelle_from_reads.py -1 ../../testdata/minimal/ptilo.R1.fq.gz -2 ../../testdata/minimal/ptilo.R2.fq.gz -F animal_mt -o getorganelle_ptilo_test/ 2>&1 >cmd_stdout | grep \"^ERROR\" | head -1"

# Try and run
errmsg=$(eval $cmd)

# If error, fix:
while [[ ! -z $errmsg ]]
do
   # Database(s) are missing
   missingdbregex="^ERROR: default (.*) database not added yet"
   if [[ $errmsg =~ $missingdbregex ]]
   then
      missingdb="${BASH_REMATCH[1]}"
      
      # Ideal command, but fails at the moment
      # get_organelle_config.py -a "$missingdb"
      
      # Workaround
      wget -O master.zip https://github.com/Kinggerm/GetOrganelleDB/archive/master.zip
      unzip -o master.zip
      get_organelle_config.py -a "$missingdb" --use-local ./GetOrganelleDB-master/0.0.1/
      
   else
      echo "$errmsg" 1>&2
      exit 1
   fi
   
   errmsg=$(eval $cmd)
done


### Notes

## Parameters
#see https://github.com/Kinggerm/GetOrganelle#starting-from-reads, especially seeds
# threads should read values passed through
# Ollie uses -s, --genes, --reduce-reads-for-coverage inf --max-reads inf

# Need to specify seed numbers for replicability
