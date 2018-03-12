#!/bin/bash
set -e

program_name="hitime-docker.sh"
input_filepath=""
output_filepath=""

# Help message for using the program.
function show_help {
cat << UsageMessage

${program_name}: detect twin ion signals in Mass Spectrometry data 

Usage:
    ${program_name} [-h] -i input.mzML -o output.mzML 

XXX description here

UsageMessage
}

# echo an error message $1 and exit with status $2
function exit_with_error {
    printf "${program_name}: ERROR: $1\n"
    exit $2
}

absolute_filepath () {
  case "$1" in
    /*) echo "$1";;
    *) echo "$PWD/$1";;
  esac
}

# Parse the command line arguments and set the global variables language and new_project_name
function parse_args {
    local OPTIND opt

    while getopts "hi:o:" opt; do
        case "${opt}" in
            h)
                show_help
                exit 0
                ;;
            i)  input_filepath="${OPTARG}"
                ;;
            o)  output_filepath="${OPTARG}"
                ;;
        esac
    done

    shift $((OPTIND-1))

    [ "$1" = "--" ] && shift

    if [[ -z ${input_filepath} ]]; then
        exit_with_error "missing command line argument: -i input.mzML, use -h for help" 2
    fi

    if [[ -z ${output_filepath} ]]; then
        exit_with_error "missing command line argument: -o output.mzML, use -h for help" 2
    fi
}

parse_args $@

input_dir=$( absolute_filepath $( dirname "${input_filepath}" ))
output_dir=$(absolute_filepath $( dirname "${output_filepath}" ))
input_filename=$( basename "${input_filepath}" )
output_filename=$( basename "${output_filepath}" )

exec docker run --rm -v "${input_dir}:/input/" -v "${output_dir}:/output/" bjpop/hitime score -i "/input/${input_filename}" -o "/output/${output_filename}"
