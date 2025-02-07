# Description: An simple script for running holmes on a single input json file.
set -e

usage() {
    echo "Usage: $0 [-j <input_json>] [-o <output_base>] [-p <holmes_path>] [-c <vep_config>] [-d <db_config>] [-t <thread_num>] [-v <vep_cache>] [-g <grch version>]" 1>&2
    exit 1
}

while getopts ":j:o:p:c:t:v:d:g:" opt; do
    case ${opt} in
        j )
            INPUT_JSON=$OPTARG
            ;;
        o )
            OUTPUT_BASE=$OPTARG
            ;;
        p )
            HOLMES=$OPTARG
            ;;
        c )
            VEP_CONFIG=$OPTARG
            ;;
        d )
            DB_CONFIG=$OPTARG
            ;;
        t )
            THREAD_NUM=$OPTARG
            ;;
        v )
            VEP_CACHE=$OPTARG
            ;;
        g )
            GRCH_VER=$OPTARG
            ;;
        \? )
            usage
            ;;
        : )
            echo "Invalid option: $OPTARG requires an argument" 1>&2
            usage
            ;;
    esac
done
shift $((OPTIND -1))

# if any of the required arguments are missing, print usage
if [ -z "$INPUT_JSON" ] || [ -z "$OUTPUT_BASE" ] || [ -z "$HOLMES" ] || [ -z "$VEP_CONFIG" ] || [ -z "$DB_CONFIG" ]; then
    usage
fi

# if thread_num is not provided, set to 8
if [ -z "$THREAD_NUM" ]; then
    echo "thread_num is not provided, set to 8"
    THREAD_NUM=8
fi

# if grch is not provided, set to 38
if [ -z "$GRCH_VER" ]; then
    echo "grch version is not provided, set to 38"
    GRCH_VER="38"
fi

# set assembly opt by $GRCH_VER
if [ "$GRCH_VER" = "38" ]; then
    ASSEMBLY_OPT=""
elif [ "$GRCH_VER" = "37" ]; then
    ASSEMBLY_OPT="--grch37"
else
    echo "grch version: $GRCH_VER is not valid, set to 38"
    GRCH_VER="38"
    ASSEMBLY_OPT=""
fi

OTHER_ARGS=${OTHER_ARGS:-""}

OUTPUT=$OUTPUT_BASE/output
VEP_OUTPUT=$OUTPUT_BASE/vep_output
LOG=$OUTPUT_BASE/log

# if vep_cache is not provided, don't include in options
if [ -z "$VEP_CACHE" ]; then
    echo "vep_cache is not provided, don't include in options"
    OPTIONS="-t $THREAD_NUM -u -d -f -j $INPUT_JSON -v $VEP_OUTPUT -o $OUTPUT --vep_config $VEP_CONFIG --db_config $DB_CONFIG $ASSEMBLY_OPT $OTHER_ARGS"
else
    OPTIONS="-t $THREAD_NUM -u -d -f -j $INPUT_JSON -v $VEP_OUTPUT -o $OUTPUT --vep_config $VEP_CONFIG --db_config $DB_CONFIG --vep_cache $VEP_CACHE $ASSEMBLY_OPT $OTHER_ARGS"
fi

mkdir -p $OUTPUT
mkdir -p $VEP_OUTPUT
mkdir -p $LOG

# run holmes
echo "[run_one.sh] Running holmes with options: $OPTIONS"
$HOLMES $OPTIONS | tee $LOG/`basename $INPUT_JSON`.log
