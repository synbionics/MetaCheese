#!/bin/bash
set -euo pipefail

###############################################################################
# Function that generates pipeline scripts from templates and config
###############################################################################
function generate_scripts() {
    local CONFIG_FILE="$1"
    local PROJECT_ROOT="$2"
    
    # Clean pipeline folder (remove old scripts)
    PIPELINE_DIR="$PROJECT_ROOT/scripts/pipeline"
    echo "Cleaning pipeline folder: $PIPELINE_DIR"
    rm -rf "$PIPELINE_DIR"/*
    mkdir -p "$PIPELINE_DIR"

    # FORMAT OF STEP DEFINITIONS:
    # [stepXX]="<TEMPLATE_FILE> <OUTPUT_FILE> <PREFIX> <ph1=YAMLkey1> <ph2=YAMLkey2> ..."
    #
    # - TEMPLATE_FILE  = name of the template script (inside scripts/templates/)
    # - OUTPUT_FILE    = path where the generated script will be written
    # - PREFIX         = prefix used in placeholders (usually the step number, e.g. 03 or 00-01)
    # - ph=YAMLkey     = mapping between:
    #                     • ph  = placeholder name used in the template (@<PREFIX>_<ph>@)
    #                     • YAMLkey = key name inside config.yml (under section stepXX)
    #
    # HOW TO CHANGE NAME OR VALUE:
    # - If you only want to change the value ? edit config.yml
    # - If you want to rename the placeholder in the template ? also change the left part (ph) in STEPS
    # - If you want to rename the key in config.yml ? also change the right part (YAMLkey) in STEPS
    #
    # Example of renaming:
    #   Template:   @03_vartemplate@
    #   STEPS:      vartemplate=newvarconfig
    #   Config.yml: step03.newvarconfig: /new/path
    # ->Result: @03_vartemplate@ becomes /new/path
    declare -A STEPS=(
        [step00_01]="template_00-01.extract_and_adapter.sh      scripts/pipeline/00-01.extract_and_adapter.sh   00-01       dir1=dir1 dir2=dir2 dir3=dir3 dir4=dir4 dir5=dir5 par1=par1 par2=par2 par3=par3 par4=par4 par5=par5 par6=par6"
        [step03]="template_03.bowtie-remove_host.sh             scripts/pipeline/03.bowtie-remove_host.sh       03          var1=dir1 var2=dir2 var3=dir3"
        [step04]="template_04.metaphlan.sh                      scripts/pipeline/04.metaphlan.sh                04          var2=dir2 var3=dir3 Rscript=rscript"
        [step04b_1]="template_04b-1.humann_job_array.sh         scripts/pipeline/04b-1.humann_job_array.sh      04b_1       var1=dir1 var2=dir2 var3=dir3"
        [step04b_last]="template_04b-last.humann_join_tables.sh scripts/pipeline/04b-last.humann_join_tables.sh 04b_last    dir1=dir1 dir2=dir2"
        [step05]="template_05.spades.assembly.sh                scripts/pipeline/05.spades.assembly.sh          05          var1=var1 var2=var2 par1=par1 par2=par2"
        [step06]="template_06.contig_filter.sh                  scripts/pipeline/06.contig_filter.sh            06          var1=var1 var2=var2 var3=var3 filter1=filter1"
        [step07]="template_07.bowtie-MAG_database.sh            scripts/pipeline/07.bowtie-MAG_database.sh      07          var1=dir1 var2=dir2 par1=par1"
        [step08_09]="template_08-09.mapp_coverage.sh            scripts/pipeline/08-09.mapp_coverage.sh         08-09       var1=dir1 var2=dir2 var3=dir3"
        [step10_11]="template_10-11-11b.metabat.sh              scripts/pipeline/10-11-11b.metabat.sh           10-11       var1=dir1 var2=dir2 var3=dir3 var4=dir4 var5=var5 par1=par1 par2=par2"
        [step12_13]="template_12-13.checkm.sh                   scripts/pipeline/12-13.checkm.sh                12-13       var1=var1 var2=var2 var3=var3 par1=par1 par2=par2 par3=par3"
        [step14]="template_14.filter-metadata.sh                scripts/pipeline/14.filter-metadata.sh          14          var1=var1 var2=var2 var3=var3 var4=var4 var5=var5 var6=var6 var7=var7 var8=var8"
        [step15]="template_15.tormes_MAG.sh                     scripts/pipeline/15.tormes_MAG.sh               15          var1=var1 var2=var2 par1=par1 par2=par2"
    )

    # Separator after prefix (usually _)
    declare -A SEPS=(
        [step00_01]="_"
        [step03]="_"
        [step04]="_"
        [step04b_1]="_"
        [step04b_last]="_"
        [step05]="_"
        [step06]="_"
        [step07]="_"
        [step08_09]="_"
        [step10_11]="_"
        [step12_13]="_"
        [step14]="_"
        [step15]="_"
    )

    # Loop through all steps
    for STEP_KEY in "${!STEPS[@]}"; do
        read -r -a arr <<< "${STEPS[$STEP_KEY]}"
        TEMPLATE_FILE="${arr[0]}"
        OUTPUT_FILE="${arr[1]}"
        PREFIX="${arr[2]}"
        map_pairs=("${arr[@]:3}")
        SEP="${SEPS[$STEP_KEY]}"

        TEMPLATE_PATH="$PROJECT_ROOT/scripts/templates/$TEMPLATE_FILE"
        OUTPUT_PATH="$PROJECT_ROOT/$OUTPUT_FILE"

        [ -f "$TEMPLATE_PATH" ] || { echo "Missing template: $TEMPLATE_PATH" >&2; continue; }
        [ -f "$CONFIG_FILE" ] || { echo "Missing config: $CONFIG_FILE" >&2; exit 1; }

        mkdir -p "$(dirname "$OUTPUT_PATH")"
        cp "$TEMPLATE_PATH" "$OUTPUT_PATH"

        # Replace placeholders @xxx@ with values taken from config.yml
        for pair in "${map_pairs[@]}"; do
            ph="${pair%%=*}"
            yaml_key="${pair#*=}"
            value=$(yq eval -r ".${STEP_KEY}.${yaml_key}" "$CONFIG_FILE")
            sed -i "s|@$PREFIX$SEP$ph@|$value|g" "$OUTPUT_PATH"
        done

        chmod +x "$OUTPUT_PATH"
    done
}

###############################################################################
# FIXED ORDER OF STEPS TO EXECUTE (modify order here if needed)
###############################################################################
STEPS_ORDERED=(
    00-01.extract_and_adapter.sh
    03.bowtie-remove_host.sh
    04.metaphlan.sh
    04b-1.humann_job_array.sh
    04b-last.humann_join_tables.sh
    05.spades.assembly.sh
    06.contig_filter.sh
    07.bowtie-MAG_database.sh
    08-09.mapp_coverage.sh
    10-11-11b.metabat.sh
    12-13.checkm.sh
    14.filter-metadata.sh
    15.tormes_MAG.sh
)

###############################################################################
# MAIN SCRIPT: Interactive menu and execution mode handling
###############################################################################

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_TEMPLATE="$PROJECT_ROOT/config/Config.yml"

# Check that yq is available
command -v yq >/dev/null 2>&1 || { echo "Error: 'yq' not found in PATH." >&2; exit 1; }

# Initial menu
echo "#########################################"
echo "           PIPELINE EXECUTION            "
echo "#########################################"
echo "Choose a mode:"
echo "1) New execution"
echo "2) Continue existing execution"
read -p "Enter 1 or 2: " scelta
while [[ "$scelta" != "1" && "$scelta" != "2" ]]; do
    read -p "Invalid choice. Enter 1 or 2: " scelta
done

read -p "Do you want to run stepwise(1) or to run the entire pipeline(2)?"
read -p "Enter 1 or 2: " risp 
while [[ "$risp" != "1" && "$risp" != "2" ]]; do
    read -p "Invalid choice. Enter 1 or 2: " risp
done
############################
# OPTION 1: NEW EXECUTION
############################
if [[ "$scelta" == "1" ]]; then
    echo "#########################################"
    echo "          NEW PIPELINE EXECUTION         "
    echo "#########################################"

    read -p "Enter input sample name: " input_sample
    [[ -z "$input_sample" ]] && { echo "Error: sample name cannot be empty!" >&2; exit 1; }

    data=$(date +%Y%m%d)
    read -p "Enter descriptive code: " codice
    [[ -z "$codice" ]] && { echo "Error: descriptive code cannot be empty!" >&2; exit 1; }

    output_folder="${data}_${codice}"
    mkdir -p "$PROJECT_ROOT/output/$output_folder"
    CONFIG_FILE="$PROJECT_ROOT/output/$output_folder/config.yml"

    # Generate config file for this run
    sed -e "s|@output_folder@|$output_folder|g" -e "s|@input_sample@|$input_sample|g" "$CONFIG_TEMPLATE" > "$CONFIG_FILE"

    # Generate pipeline scripts from templates
    generate_scripts "$CONFIG_FILE" "$PROJECT_ROOT"

    echo "Pipeline successfully generated. Output folder: $output_folder"

    cd "$PROJECT_ROOT/scripts/pipeline"

    # Execute ALL scripts in STEPS_ORDERED
    for step_script in "${STEPS_ORDERED[@]}"; do
        if [ -f "$step_script" ]; then
            echo "***************************"
            echo " Executing $step_script"
            echo "***************************"
            bash "$step_script" || { echo "Error in $step_script, stopping." >&2; exit 1; }
            echo "***************************"
            echo " Completed $step_script"
            echo "***************************"
            #read -p "Do you want to continue with the next step? (y/n): " risp
            [[ "$risp" =~ ^[Nn]$ ]] && { echo "Pipeline stopped on user request."; exit 0; }
        else
            echo "Warning: $step_script not found! Skipping."
        fi
    done

    echo "Pipeline successfully completed."

############################
# OPTION 2: CONTINUE EXECUTION
############################
elif [[ "$scelta" == "2" ]]; then
    read -p "Enter existing output_folder name: " output_folder
    OUTPUT_DIR="$PROJECT_ROOT/output/$output_folder"
    CONFIG_FILE="$OUTPUT_DIR/config.yml"

    [[ ! -d "$OUTPUT_DIR" ]] && { echo "Error: folder $OUTPUT_DIR does not exist!" >&2; exit 1; }
    [[ ! -f "$CONFIG_FILE" ]] && { echo "Error: config.yml not found in $OUTPUT_DIR" >&2; exit 1; }

    generate_scripts "$CONFIG_FILE" "$PROJECT_ROOT"
    echo "Scripts regenerated from $CONFIG_FILE"

    PIPELINE_DIR="$PROJECT_ROOT/scripts/pipeline"
    echo "Available scripts:"
    ls "$PIPELINE_DIR"/*.sh | xargs -n1 basename
    read -p "Enter part of the script name to resume from: " search_term

    # Find index of the script to resume from
    start_index=-1
    for i in "${!STEPS_ORDERED[@]}"; do
        if [[ "${STEPS_ORDERED[$i]}" == *"$search_term"* ]]; then
            start_index=$i
            break
        fi
    done

    if [[ $start_index -eq -1 ]]; then
        echo "No script found with '$search_term'"
        exit 1
    fi

    cd "$PIPELINE_DIR"
    # Resume execution from the correct step, maintaining the order
    for (( i=$start_index; i<${#STEPS_ORDERED[@]}; i++ )); do
        step_script="${STEPS_ORDERED[$i]}"
        if [ -f "$step_script" ]; then
            echo "***************************"
            echo " Executing $step_script"
            echo "***************************"
            bash "$step_script" || { echo "Error in $step_script" >&2; exit 1; }
            echo "***************************"
            echo " Completed $step_script"
            echo "***************************"
            #read -p "Do you want to continue with the next step? (y/n): " risp
            [[ "$risp" =~ ^[Nn]$ ]] && { echo "Pipeline stopped on user request."; exit 0; }
        else
            echo "Warning: $step_script not found! Skipping."
        fi
    done
    echo "Pipeline successfully completed."
fi