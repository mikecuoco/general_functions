#!/usr/bin/awk -f

# Note: awk proceeds line-by-line, so it goes through all of this code for the first
#       line and then goes through it all again for the second line, etc

# Print header lines without further processing
# We use a regex to execute the code in the brackets only for header lines which we
# know will begin with the @ symbol
/^@/ {
    # send the line to stdout
    print;
    # skip the rest of the code in this script when processing this line
    # this will force awk to move on to the next line
    next;
}

# This code will be executed for all other lines:
{
    # lots of stuff happening in this line:
    # First, we create an array and store it in the variable 'counts'
    # Arrays in awk are associative. So they're really just key-value pairs $1 refers
    # to the value in the first column of the line, and by doing counts[$1], we're
    # just accessing the value of the array with key $1
    # Note that awk array values automatically initialize to 0
    # The ++ part returns the value and then increments it, so if the read ID in $1
    # hasn't been seen before, then counts[$1]++ will be 0 upon the first time we see
    # it. It will be 1 upon the next time we see it, and 2 the following time, etc.
    # Once we've acquired the number of times we've seen this read ID before, we store
    # it in the variable 'count'
    count = counts[$1]++
    # If count is true (ie not 0)...
    if (count) {
        b = 0;
    } else {
        if (b) {
            print a
        }
        # If this is the first time we've seen this read ID, we store the entirety of
        # the line in a variable called 'a', in case we need it later.
        a = $0;
        b = 1;
    }
}

END{
    if (b) {
        print a
    }
}