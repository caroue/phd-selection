#!usr/bin/bash

while getopts u:a:f: flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        o) output=${OPTARG};;
    esac
done

awk '/### Adaptive branch site random effects likelihood test/,0' $input |\
grep -i "^\*" |\
awk -F ' ' '{print $2}' |\
sed -i 's/,//g'
> $output
