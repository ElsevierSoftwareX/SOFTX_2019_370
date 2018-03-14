#!/bin/bash
set -e

program_name="hitime-docker.sh"
input_filepath=""
output_filepath=""
# this should be one of "score" or "max"
mode=""
remaining_args=""

# Help message for using the program.
function show_help {
cat << UsageMessage

${program_name}: detect twin ion signals in Mass Spectrometry data 

Usage:
    ${program_name} [-h] (score|max) -i input.mzML -o output.mzML [-- optional mode-specific arguments]

This is a wrapper for the HiTIME-CPP docker container.

You must choose a mode, either score or max.

Both modes read and write mzML files.

The score mode re-scales the peaks in the input file based on their similarity to
an ideal twin ion peak.

The max mode computes a local maxima for peaks within a defined window. 

Example usage:

   ./hitime-docker.sh score -i data/testing.mzML -o data/scored.mzML -- -j 2

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

    if [ "$#" -lt 1 ]; then
        show_help
        exit 2 
    fi

    if [ "$1" = 'score' ]; then
        shift 1
        mode='score'
    elif [ "$1" = 'max' ]; then
        shift 1
        mode='max'
    fi 

    local OPTIND opt

    while getopts "hi:o:" opt; do
        case "${opt}" in
            h)
                # Print mode specific help if mode is defined
                if [[ -n ${mode} ]]; then
                    exec docker run --rm bjpop/hitime ${mode} -h
                fi
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

    remaining_args=$@

    if [[ -z ${mode} ]]; then
        exit_with_error "mode not specified, must be one of 'score' or 'max', use -h for help" 2
    fi

    if [[ -z ${input_filepath} ]]; then
        exit_with_error "input file not specified: -i <filename>, use -h for help" 2
    fi

    if [[ -z ${output_filepath} ]]; then
        exit_with_error "output file not specified: -o <filename>, use -h for help" 2
    fi

}

parse_args $@

input_dir=$( absolute_filepath $( dirname "${input_filepath}" ))
output_dir=$(absolute_filepath $( dirname "${output_filepath}" ))
input_filename=$( basename "${input_filepath}" )
output_filename=$( basename "${output_filepath}" )

exec docker run --rm -v "${input_dir}:/input/" -v "${output_dir}:/output/" bjpop/hitime ${mode} -i "/input/${input_filename}" -o "/output/${output_filename}" ${remaining_args}
