#!/bin/bash
# Run all test combinations

# Run centralized tests

main() {
    tests/centralized-parallel-timing setup-cent-par
    sleep 5

    tests/centralized-serial-timing setup-cent-ser
    sleep 5

    for i in $(seq 1 9); do
        echo "Using $i timing iterations."
        n_timing_iterations=$i

        filename="coop$i.dat"
        set_setup_params "setup-coop-par"
        tests/cooperative-parallel-timing setup-coop-par
        sleep 5

        set_setup_params "setup-coop-ser"
        tests/cooperative-serial-timing setup-coop-ser
        sleep 5

        filename="ncoop$i.dat"
        set_setup_params "setup-ncoop-par"
        tests/noncoop-parallel-timing setup-ncoop-par
        sleep 5

        set_setup_params "setup-ncoop-ser"
        tests/noncoop-serial-timing setup-ncoop-ser
        sleep 5

    done

}


set_setup_params() {
    mv "$1" setup.tmp
    gawk 'BEGIN { found=0; } /n-timing-iterations/ { found=1; print; next; } { if(found) { found=0; print n; next; } } 1' n="$n_timing_iterations" setup.tmp > "$1"
    mv "$1" setup.tmp

    gawk 'BEGIN { found=0; } /output-filename/ { found=1; print; next; } { if(found) { found=0; print name; next; } } 1' name="$filename" setup.tmp > "$1"
    rm setup.tmp
}

main

