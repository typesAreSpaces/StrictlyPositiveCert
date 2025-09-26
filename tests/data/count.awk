#!/opt/homebrew/bin/awk -f

BEGIN {
  count = 0
}

/^\[\[\], \[\], \[\]\]/ {
    count = count + 1
}

END {
    print "Number of successful cases " count
}
