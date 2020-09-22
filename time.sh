#!/bin/bash
START=$(date +%N)
# do something
# start your script work here
./client 5
# your logic ends here
END=$(date +%N)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds"